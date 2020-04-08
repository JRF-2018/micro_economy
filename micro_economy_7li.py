#!/usr/bin/python3
__version__ = '0.0.10' # Time-stamp: <2020-03-10T05:24:33Z>
## Language: Japanese/UTF-8

"""A simple simulation of micro economics."""

##
## License:
##
##   Public Domain
##   (Since this small code is close to be mathematically trivial.)
##
## Author:
##
##   JRF
##   http://jrf.cocolog-nifty.com/software/
##   (The page is written in Japanese.)
##

import numpy as np
import random
import scipy.optimize
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import chi2
import traceback
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--trials", default=20, type=int)
parser.add_argument("--opt-try", default="none", choices=["none", "prev", "random", "bh0", "de", "shgo", "da", "da0"])
parser.add_argument("--score-labor-mag", default=0.00005, type=float)
parser.add_argument("--score-mean-mag", default=0.0, type=float)
parser.add_argument("--score-variance-mag", default=10.0, type=float)
parser.add_argument("--score-variance-mag-2", default=1.0, type=float)
parser.add_argument("--score-transform", default="none", choices=["none", "log", "log1.1", "exp", "hybrid", "hybrid0.2", "hybrid0.1"])
parser.add_argument("--score-mean", default="max", choices=["max", "mean", "each"])
parser.add_argument("--expected-companies", default=5.0, type=float)
parser.add_argument("--score-companies-mag", default=1, type=float)
parser.add_argument("--score-profit-rate-mag", default=0.01, type=float)
parser.add_argument("--score-social-debt-mag", default=1000.0, type=float)
parser.add_argument("--score-zero-worker-penalty", default=100.0, type=float)
parser.add_argument("--init-lpp-mag", default=1.0, type=float)
parser.add_argument("--fix-standard-profit-rate", default=False, action="store_true")
parser.add_argument("--expected-standard-profit-rate", default=0.06, type=float)
parser.add_argument("--necessaries-elasticity", default=-0.6, type=float)
parser.add_argument("--surplus-vs-necessaries-ratio", default=0.5, type=float)
parser.add_argument("--wages-elastic-power", default=1.0, type=float)
parser.add_argument("--wages-elastic-power-2", default=1.0, type=float)
parser.add_argument("--ingredients-elastic-power", default=1.0, type=float)
parser.add_argument("--ingredients-elastic-power-2", default=1.0, type=float)
parser.add_argument("--ingredients-limit", default=0.9, type=float)
parser.add_argument("--score-ingredients-level-mag", default=0.01, type=float)
parser.add_argument("--score-wrong-ingredients-level-mag", default=0.01, type=float)
parser.add_argument("--power-necessaries-ingredients", default=0.0, type=float)
parser.add_argument("--update-workers-level", default=1.0, type=float)
parser.add_argument("--save-history", default=None, type=str)

ARGS = parser.parse_args()


def save_history (path, history):
    import csv
    with open(path, 'w') as f:
        keys = list(history.keys())
        epochs = len(history[keys[0]])
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        r = [dict([[k, history[k][i]] for k in keys])
             for i in range(epochs)]
        writer.writerows(r)


class Labor:
    def __init__ (self):
        self.working = True
        self.savings = None

class Capitalist:
    def __init__ (self):
        self.id = None
        self.working = True
        self.savings = None
        self.debt = None
        self.makeCompany = self._makeCompany_r

    def _makeCompany (self, oldCompany, curPrice, igLimit, r1, r2):
        c = Company()
        o = oldCompany
        c.ownerId = self.id
        c.fixedAssetsIngredients = o.fixedAssetsIngredients * 1.05
        c.fixedAssetsPrice = c.fixedAssetsIngredients * curPrice.ingredients
        c.fixedAssetsLife = 10
        c.fixedAssetsLifeRest = c.fixedAssetsLife
        c.superiority = (4 - 2) * r1 + 2
        c.superiorityDecay = ((0.5 ** 0.1) - (0.25 ** 0.1)) * r2 + (0.25 ** 0.1)
        c.laborsPerProduct = ((0.99 - 0.95) * r1 + 0.95) * o.laborsPerProduct
        r = ((1.18 - 0.98) * r2 + 0.98) * o.ingredientsPerProduct
        if igLimit is not None and r > igLimit:
            r = igLimit
        c.ingredientsPerProduct = r
        return c

    def _makeCompany_r (self, oldCompany, curPrice, igLimit):
        return self._makeCompany(oldCompany, curPrice, igLimit, \
                                     random.random(), random.random())
    def _makeCompany_0_0 (self, oldCompany, curPrice, igLimit):
        return self._makeCompany(oldCompany, curPrice, igLimit, 0, 0)
    def _makeCompany_0_1 (self, oldCompany, curPrice, igLimit):
        return self._makeCompany(oldCompany, curPrice, igLimit, 0, 1)
    def _makeCompany_1_0 (self, oldCompany, curPrice, igLimit):
        return self._makeCompany(oldCompany, curPrice, igLimit, 1, 0)
    def _makeCompany_1_1 (self, oldCompany, curPrice, igLimit):
        return self._makeCompany(oldCompany, curPrice, igLimit, 1, 1)

class Company:
    def __init__ (self):
        self.ownerId = None
        self.fixedAssetsPrice = None
        self.fixedAssetsIngredients = None
        self.fixedAssetsLife = None
        self.fixedAssetsLifeRest = None
        self.superiority = None
        self.superiorityDecay = None
        self.laborsPerProduct = None
        self.ingredientsPerProduct = None

class Price:
    def __init__ (self, priceArray):
        self.wages = priceArray[0]
        self.necessaries = priceArray[1]
        self.luxuries = priceArray[2]
        self.ingredients = priceArray[3]
        self.standardProfitRate = priceArray[4]
        self.ingredientsLevelForLuxuries = priceArray[5]

    def toArray (self):
        return [self.wages, self.necessaries, self.luxuries,
                self.ingredients, self.standardProfitRate,
                self.ingredientsLevelForLuxuries]

def randomList (l):
    return random.sample(l, len(l))

def rand1to100 ():
    while True:
        x = random.uniform(1.0, 100.0)
#        x = random.uniform(1.0, 1000.0)
        if x != 1.0 and x != 1000:
            return x

def _score_transform_log(x):
    return np.log(1 + x)

def _score_transform_log_1_1(x):
    a = 1.1
    c = 1.1
    b = - a * np.log(c)
    return b + a * np.log(x + c)

def _score_transform_exp(x):
    b = 1
    c = -5  # x = -1 のときの y の値。
    a = np.log(1 - (c/b))
    return - b * np.exp(- a * x) + b

def _score_transform_hybrid(x):
    if x >= 0:
        return 0.2 * np.log(x + 0.1) - 0.2 * np.log(0.1)
    else:
        return - np.exp(- 2 * x) + 1

def _score_transform_hybrid_0_2(x):
    if x >= 0:
        return 0.4 * np.log(x + 0.2) - 0.4 * np.log(0.2)
    else:
        return - np.exp(- 2 * x) + 1

def _score_transform_hybrid_0_1(x):
    if x >= 0:
        return 0.1 * np.log(x + 0.1) - 0.1 * np.log(0.1)
    else:
        return - np.exp(- 2 * x) + 1

score_transform = {
    'none': lambda x: x,
    'log': _score_transform_log,
    'log1.1': _score_transform_log_1_1,
    'exp': _score_transform_exp,
    'hybrid': _score_transform_hybrid,
    'hybrid0.2': _score_transform_hybrid_0_2,
    'hybrid0.1': _score_transform_hybrid_0_1,
}[ARGS.score_transform]


def score_expected_companies (n):
    df = ARGS.expected_companies + 2
    return chi2.pdf(n, df) / chi2.pdf(ARGS.expected_companies, df)


# abbrev.: wL == workersLevel, iL == ingredientsLevel
class Economy:
    def __init__ (self):
        random.seed()
        self.debug = False
        self.randomstate = random.getstate()

        self.ExpectedStandardProfitRate = ARGS.expected_standard_profit_rate
        self.NecessariesElasticity = ARGS.necessaries_elasticity
        self.SurplusVsNecessariesRatio = ARGS.surplus_vs_necessaries_ratio
        self.WagesElasticPower = ARGS.wages_elastic_power
        self.WagesElasticPower2 = ARGS.wages_elastic_power_2
        self.IngredientsElasticPower = ARGS.ingredients_elastic_power
        self.IngredientsElasticPower2 = ARGS.ingredients_elastic_power_2
        self.IngredientsLimit = ARGS.ingredients_limit

        self.cumulativeSocialDebt = 0
        self.capitalists = []
        self.labors = []
        self.necessariesCompanies = []
        self.luxuriesCompanies = []
        self.ingredientsCompanies = []

        self.lastCompanyOfNecessaries = None
        self.lastCompanyOfLuxuries = None
        self.lastCompanyOfIngredients = None

        self.profitHistoryOfNecessaries = []
        self.profitHistoryOfLuxuries = []
        self.profitHistoryOfIngredients = []
        self.prevPrice = None
        self.prevDemandOfNecessaries = None
        self.prevSocialDebt = 0

        self.priceHistory = []

        self.workersLevel = 1.0
        self.lppSave = None
        self.ippSave = None

        self.history = {}
        history_labels = [
            'score',
            'price_W', 'price_N', 'price_L', 'price_I',
            'standard_profit_rate',
            'num_companies_N', 'num_companies_L', 'num_companies_I',
            'social_debt', 'savings', 'savings_increase',
            'workers_savings_increase',
            'supply_of_labors', 'demand_vs_supply_of_labors',
            'max_profit_rate_N', 'max_profit_rate_L',
            'max_profit_rate_I', 'max_profit_rate_mean',
            'max_profit_rate_all',
            'mean_profit_rate_N', 'mean_profit_rate_L',
            'mean_profit_rate_I', 'mean_profit_rate_all',
            'product_N', 'product_L', 'product_I',
            'savings_mean', 'savings_variance',
            'ingredients_level_for_luxuries', 'workers_level'
            ]

        for l in history_labels:
            self.history[l] = []

        if ARGS.opt_try != "none":
            self.history['score_1'] = []
            self.history['score_2'] = []


    def _calcProfitOfCompany (self, company, price, curPrice, num):
        c = company
        cost = (c.fixedAssetsPrice / c.fixedAssetsLife) \
            + (curPrice.wages * c.laborsPerProduct
               + curPrice.ingredients * c.ingredientsPerProduct) * num
        return (price * num - cost) / cost

    def _calcProducts (self, curPrice, price, demand, companies):
        curS = [0 for c in companies]
        delta = self.WagesElasticPower2 \
            * np.tanh(self.WagesElasticPower \
                          * (curPrice.wages - self.prevPrice.wages)
                      / self.prevPrice.wages)
        theta = self.IngredientsElasticPower2 \
            * np.tanh(self.IngredientsElasticPower \
                          * (curPrice.ingredients
                             - self.prevPrice.ingredients) \
                          / self.prevPrice.ingredients)
        for i, c in enumerate(companies):
            if (price - curPrice.wages * c.laborsPerProduct \
                    - curPrice.ingredients * c.ingredientsPerProduct) > 0:
                curS[i] = c.superiority \
                    * (c.laborsPerProduct ** (- delta)) \
                    * (c.ingredientsPerProduct ** (- theta))
        s = sum(curS)
        if s != 0:
            for i in range(len(curS)):
                curS[i] = (curS[i] / s) * demand
        else:
            for i, c in enumerate(companies):
                curS[i] = c.superiority * (c.laborsPerProduct ** (- delta))
            s = sum(curS)
            for i in range(len(companies)):
                curS[i] = (curS[i] / s) * demand
        return curS

    def _scoreOfPOI (self, k, curPrice, price, initDemand, companies):
        if k < 0:
#            return float("inf") ## なぜかエラーになる。
            return initDemand ** 2
        u = k + initDemand
        curS = self._calcProducts(curPrice, price, u, companies)
        s = initDemand
        for i, c in enumerate(companies):
            s += curS[i] * c.ingredientsPerProduct
        return (s - u) ** 2
    
    def _calcProductsOfIngredients (self, curPrice, price, initDemand,
                                    companies):
        maxI = max(map(lambda x: x.ingredientsPerProduct, companies))
        if initDemand <= 0:
            print("ERROR: initDemand must be > 0.")
        res = minimize_scalar(
            lambda k: self._scoreOfPOI(k, curPrice, price, 
                                       initDemand, companies), \
                bracket=(0, (maxI / (1 - maxI)) * initDemand), method='brent')
        sc = self._scoreOfPOI(res.x, curPrice, price, initDemand, companies)
        if sc >= 0.1:
            curS = self._calcProducts(curPrice, price, res.x + initDemand, 
                                      companies)
            print(("Warning: Bad score about ingredients: score = {0}, "
                   + "k = {1}, d0 = {2}, iterated {3} times.")
                  .format(sc, res.x, initDemand, res.nit))
        return self._calcProducts(curPrice, price, res.x + initDemand,
                                  companies)

    def _cp_debt_key (self, cp):
        if cp.debt >= 0.1:
            return 0
        else:
            return 1

    def _calcNextInitialValue(self, clip=True):
        l = self.prevPrice.toArray()
        if ARGS.fix_standard_profit_rate:
            p = self.ExpectedStandardProfitRate
        else:
            p = np.clip(
                np.mean(list(map(lambda l: l[-1],
                                 [self.profitHistoryOfNecessaries,
                                  self.profitHistoryOfLuxuries,
                                  self.profitHistoryOfIngredients]))),
                0.0, 1.0)
        if clip:
            for i in range(4):
                l[i] = np.clip(l[i], 1, 100)
            l[5] = np.clip(l[i], 1, 100)

        return l[0:4] + [p] + l[5:]

    def _setTemporalLcIpp(self, iL, wL):
        if self.ippSave is not None:
            raise ValueError("originalLcIpp has already been set.")
        nL = iL ** ARGS.power_necessaries_ingredients
        self.lppSave = {}
        self.ippSave = {}
        l1, l2 = [], []
        for c in self.necessariesCompanies:
            l1.append(c.laborsPerProduct)
            l2.append(c.ingredientsPerProduct)
            c.laborsPerProduct *= wL
            c.ingredientsPerProduct *= nL
        self.lppSave['necessaries'] = l1
        self.ippSave['necessaries'] = l2

        l1, l2 = [], []
        for c in self.luxuriesCompanies:
            l1.append(c.laborsPerProduct)
            l2.append(c.ingredientsPerProduct)
            c.laborsPerProduct *= wL
            c.ingredientsPerProduct *= iL
        self.lppSave['luxuries'] = l1
        self.ippSave['luxuries'] = l2

        l1, l2 = [], []
        for c in self.ingredientsCompanies:
            l1.append(c.laborsPerProduct)
            l2.append(c.ingredientsPerProduct)
            c.laborsPerProduct *= wL
        self.lppSave['ingredients'] = l1
        self.ippSave['ingredients'] = l2

    def _restoreLcIpp(self):
        for n, companies in zip(['necessaries', 'luxuries', 'ingredients'],
                                [self.necessariesCompanies, 
                                 self.luxuriesCompanies,
                                 self.ingredientsCompanies]):
            l1 = self.lppSave[n]
            l2 = self.ippSave[n]
            for c, x1, x2 in zip(companies, l1, l2):
                c.laborsPerProduct = x1
                c.ingredientsPerProduct = x2
        self.lppSave = None
        self.ippSave = None
            

    def _calcNextState (self, priceArray):
        random.setstate(self.randomstate)

        curPrice = Price(priceArray)

        l = priceArray
        if any(map(lambda x: not (x >= 1 and x <= 1000), l[0:4])):
            return {'score': float('inf')}

        if not (l[4] >= 0.0 and l[4] <= 1.0):
            return {'score': float('inf')}

        if not (l[5] > 0.0 and l[5] <= 100):
            return {'score': float('inf')}

        if ARGS.fix_standard_profit_rate:
            standardProfitRate = self.ExpectedStandardProfitRate
        else:
            standardProfitRate = curPrice.standardProfitRate

        self._setTemporalLcIpp(curPrice.ingredientsLevelForLuxuries,
                               self.workersLevel)


        ## 必需品需要の計算
        prevD = self.prevDemandOfNecessaries \
            / (len(self.capitalists) + len(self.labors))
        curD = (1 + ((curPrice.necessaries - self.prevPrice.necessaries)
                     / self.prevPrice.necessaries)
                * self.NecessariesElasticity) * prevD
        if curD < 1:
            curD = 1
        curDemandOfNecessaries = curD \
            * (len(self.capitalists) + len(self.labors))

        ## 労働供給の計算
        curSupplyOfLabors = None
        curS = curPrice.wages - curPrice.necessaries * curD
        curSavingsOfLabors = curS
        ## チューンアップ前: (新規貯蓄をスコアに加えてもうまくいかない…。)
        if curS <= 0:
        ## チューンアップ後:
#        if curS <= curPrice.wages / 3:
            curSupplyOfLabors = 0
        else:
            prevS = self.prevPrice.wages - self.prevPrice.necessaries * prevD
            q = ((curS - prevS) / prevS) * self.SurplusVsNecessariesRatio \
                + ((curPrice.necessaries - self.prevPrice.necessaries) \
                       / self.prevPrice.necessaries) \
                       * (1 - self.SurplusVsNecessariesRatio)
            r = 0.5 * 2 * np.arctan(2 * q) / np.pi
            working = 0
            for l in self.labors:
                if l.working:
                    working += 1
            notWorking = len(self.labors) - working
            if notWorking > working:
                curSupplyOfLabors = working + np.ceil(r * working)
            else:
                curSupplyOfLabors = len(self.labors) \
                    - (notWorking - np.floor(r * notWorking))

        ## 贅沢品需要の計算
        base = curS * 3.0
        if base < 0:
            base = 0
        curDemandOfLuxuries = 0
        for i, x in enumerate(self.labors + self.capitalists):
            if i >= len(self.labors) or i < curSupplyOfLabors:
                if x.savings + curS >= base \
                        and ((x.savings + curS - base) / 3) + base \
                        >= curPrice.luxuries:
                    curDemandOfLuxuries += 1
            else:
                if x.savings >= base \
                        and ((x.savings - base) / 3) + base \
                        >= curPrice.luxuries:
                    curDemandOfLuxuries += 1
                

        ## 必需品生産
        newC = None
        companies = self.necessariesCompanies
        prHistory = self.profitHistoryOfNecessaries
        lastC = self.lastCompanyOfNecessaries
        p = curPrice.necessaries
        u = curDemandOfNecessaries
        for i, cp in enumerate(sorted(randomList(self.capitalists),
                                      key=self._cp_debt_key)):
            ## 投資のためには debt がないのが条件
            if cp.debt >= 0.1 \
                    and not (len(companies) == 0 \
                                 and i == len(self.capitalists) - 1):
                continue
            self._restoreLcIpp()
            newC = cp.makeCompany(lastC, curPrice, None)
            companies.append(newC)
            self._setTemporalLcIpp(curPrice.ingredientsLevelForLuxuries,
                                   self.workersLevel)
            curP = self._calcProducts(curPrice, p, u, companies)
            ## その期の利潤率と3期平均利潤率の低い方が市場利子率より
            ## 高いとき投資する。
            pr = self._calcProfitOfCompany(newC, p, curPrice, curP[-1])
            pr2 = sum(prHistory[-2:] + [pr]) / 3
            if pr > pr2:
                pr = pr2
            ## 会社が 0 の場合は最後の資本家が責任を持つ。
            if  pr >= standardProfitRate \
                    or (len(companies) <= 1 \
                            and i == len(self.capitalists) - 1):
                break
            else:
                newC = None
                companies.pop()
        if newC is None:
            curP = self._calcProducts(curPrice, p, u, companies)
        curProductsOfNecessaries = curP
        newCompanyOfNecessaries = newC


        ## 贅沢品生産
        newC = None
        curP = None
        companies = self.luxuriesCompanies
        prHistory = self.profitHistoryOfLuxuries
        lastC = self.lastCompanyOfLuxuries
        p = curPrice.luxuries
        u = curDemandOfLuxuries
        for i, cp in enumerate(sorted(randomList(self.capitalists),
                                      key=self._cp_debt_key)):
            ## 投資のためには debt がないのが条件
            if cp.debt >= 0.1 \
                    and not (len(companies) == 0 \
                                 and i == len(self.capitalists) - 1):
                continue
            self._restoreLcIpp()
            newC = cp.makeCompany(lastC, curPrice, None)
            companies.append(newC)
            self._setTemporalLcIpp(curPrice.ingredientsLevelForLuxuries,
                                   self.workersLevel)
            curP = self._calcProducts(curPrice, p, u, companies)
            ## その期の利潤率と3期平均利潤率の低い方が市場利子率より
            ## 高いとき投資する。
            pr = self._calcProfitOfCompany(newC, p, curPrice, curP[-1])
            pr2 = sum(prHistory[-2:] + [pr]) / 3
            if pr > pr2:
                pr = pr2
            ## 会社が 0 の場合は最後の資本家が責任を持つ。
            if  pr >= standardProfitRate \
                    or (len(companies) <= 1 \
                            and i == len(self.capitalists) - 1):
                break
            else:
                newC = None
                companies.pop()
        if newC is None:
            curP = self._calcProducts(curPrice, p, u, companies)
        curProductsOfLuxuries = curP
        newCompanyOfLuxuries = newC
        
        ## 原料生産
        newC = None
        curP = None
        companies = self.ingredientsCompanies
        prHistory = self.profitHistoryOfIngredients
        lastC = self.lastCompanyOfIngredients
        p = curPrice.ingredients
        u = 0
        for i, c in enumerate(self.necessariesCompanies):
            u += curProductsOfNecessaries[i] * c.ingredientsPerProduct
        for i, c in enumerate(self.luxuriesCompanies):
            u += curProductsOfLuxuries[i] * c.ingredientsPerProduct
        if newCompanyOfNecessaries is not None:
            u += newCompanyOfNecessaries.fixedAssetsIngredients
        if newCompanyOfLuxuries is not None:
            u += newCompanyOfLuxuries.fixedAssetsIngredients
        for i, cp in enumerate(sorted(randomList(self.capitalists),
                                      key=self._cp_debt_key)):
            ## 投資のためには debt がないのが条件
            if cp.debt >= 0.1 \
                    and not (len(companies) == 0 \
                                 and i == len(self.capitalists) - 1):
                continue
            self._restoreLcIpp()
            newC = cp.makeCompany(lastC, curPrice, self.IngredientsLimit)
            companies.append(newC)
            self._setTemporalLcIpp(curPrice.ingredientsLevelForLuxuries,
                                   self.workersLevel)
            curP = self._calcProductsOfIngredients(\
                curPrice, p, u + newC.fixedAssetsIngredients, companies)
            ## その期の利潤率と3期平均利潤率の低い方が市場利子率より
            ## 高いとき投資する。
            pr = self._calcProfitOfCompany(newC, p, curPrice, curP[-1])
            pr2 = sum(prHistory[-2:] + [pr]) / 3
            if pr > pr2:
                pr = pr2
            ## 会社が 0 の場合は最後の資本家が責任を持つ。
            if  pr >= standardProfitRate \
                    or (len(companies) <= 1 \
                            and i == len(self.capitalists) - 1):
                break
            else:
                newC = None
                companies.pop()
        if newC is None:
            if len(companies) == 0:
                print("WHY")
            curP = self._calcProductsOfIngredients(curPrice, p, u, companies)
        curProductsOfIngredients = curP
        newCompanyOfIngredients = newC

        ## 労働需要の計算
        r = 0
        curP = curProductsOfNecessaries
        for i, c in enumerate(self.necessariesCompanies):
            r += curP[i] * c.laborsPerProduct
        curP = curProductsOfLuxuries
        for i, c in enumerate(self.luxuriesCompanies):
            r += curP[i] * c.laborsPerProduct
        curP = curProductsOfIngredients
        for i, c in enumerate(self.ingredientsCompanies):
            r += curP[i] * c.laborsPerProduct
        curDemandOfLabors = r


        ## 社会的負債 sd の計算
        sd = 0

        if curDemandOfLabors > curSupplyOfLabors:
            sd += 2 * curPrice.wages * (curDemandOfLabors - curSupplyOfLabors)
        else:
            sd += (curSupplyOfLabors - curDemandOfLabors) * curPrice.wages
        curD = curDemandOfNecessaries \
            / (len(self.capitalists) + len(self.labors))
        sd += curPrice.necessaries * curD \
            * (len(self.labors) - curSupplyOfLabors)

        ## 資本家利潤
        debt = [x.debt for x in self.capitalists]
        newC = newCompanyOfNecessaries
        if newC is not None:
            debt[newC.ownerId] += newC.fixedAssetsPrice
        newC = newCompanyOfLuxuries
        if newC is not None:
            debt[newC.ownerId] += newC.fixedAssetsPrice
        newC = newCompanyOfIngredients
        if newC is not None:
            debt[newC.ownerId] += newC.fixedAssetsPrice

        profit = [0 for i in self.capitalists]
        debtable = [0 for i in self.capitalists]
        curP = curProductsOfNecessaries
        p = curPrice.necessaries
        for i, c in enumerate(self.necessariesCompanies):
            debtable[c.ownerId] += c.fixedAssetsPrice / c.fixedAssetsLife
            profit[c.ownerId] = curP[i] * (\
                p - curPrice.wages * c.laborsPerProduct \
                    + curPrice.ingredients * c.ingredientsPerProduct)

        curP = curProductsOfLuxuries
        p = curPrice.luxuries
        for i, c in enumerate(self.luxuriesCompanies):
            debtable[c.ownerId] += c.fixedAssetsPrice / c.fixedAssetsLife
            profit[c.ownerId] = curP[i] * (\
                p - curPrice.wages * c.laborsPerProduct \
                    + curPrice.ingredients * c.ingredientsPerProduct)

        curP = curProductsOfIngredients
        p = curPrice.ingredients
        for i, c in enumerate(self.ingredientsCompanies):
            debtable[c.ownerId] += c.fixedAssetsPrice / c.fixedAssetsLife
            profit[c.ownerId] = curP[i] * (\
                p - curPrice.wages * c.laborsPerProduct \
                    + curPrice.ingredients * c.ingredientsPerProduct)

        for i, cp in enumerate(self.capitalists):
            profit[i] -= curPrice.wages
            debtable[i] += curPrice.wages
            sd += debtable[i]
            debt[i] -= debtable[i]
            if debt[i] < 0:
                sd += debt[i]
                debt[i] = 0
            if profit[i] < 0:
                sd += - profit[i]
            else:
                debt[i] -= profit[i]
                if debt[i] < 0:
                    sd += debt[i]
                    debt[i] = 0


        ## スコアの計算
        r = 0

        ## 労働供給が 0 になる(新規貯蓄が 0 以下になる)のは実現不能なの
        ## で社会的負債の計算から離れペナルティを与えることにする。
        if curSupplyOfLabors == 0:
            r += ARGS.score_zero_worker_penalty

        r += np.tanh(ARGS.score_social_debt_mag
                       * ((sd - self.prevSocialDebt) / 
                          (1000 if abs(self.prevSocialDebt) < 1000
                           else abs(self.prevSocialDebt))))

        r += np.tanh(ARGS.score_labor_mag
                     * ((curDemandOfLabors - curSupplyOfLabors) ** 2))

        r0 = []

        curP = curProductsOfNecessaries
        p = curPrice.necessaries
        r1 = []
        for i, c in enumerate(self.necessariesCompanies):
            r1.append(self._calcProfitOfCompany(c, p, curPrice, curP[i]))
        r0.append(r1)

        curP = curProductsOfLuxuries
        p = curPrice.luxuries
        r1 = []
        for i, c in enumerate(self.luxuriesCompanies):
            r1.append(self._calcProfitOfCompany(c, p, curPrice, curP[i]))
        r0.append(r1)

        curP = curProductsOfIngredients
        p = curPrice.ingredients
        r1 = []
        for i, c in enumerate(self.ingredientsCompanies):
            r1.append(self._calcProfitOfCompany(c, p, curPrice, curP[i]))
        r0.append(r1)

        mm = min(map(np.mean, r0))

        if ARGS.score_mean == 'max':
            r0 = map(lambda l: score_transform(max(l)), r0)
        elif ARGS.score_mean == 'mean':
            r0 = map(lambda l: score_transform(np.mean(l)), r0)
        else:
            r0 = [score_transform(item)
                  for sublist in r0 for item in sublist]

        r0 = list(r0)
        r += - np.tanh(ARGS.score_mean_mag * np.mean(r0)) \
            + ARGS.score_variance_mag_2 \
            * np.tanh(ARGS.score_variance_mag * np.var(r0))

        r_c = ARGS.score_companies_mag * sum(
            map(lambda l: score_expected_companies(len(l)),
                [self.necessariesCompanies,
                 self.luxuriesCompanies,
                 self.ingredientsCompanies])
            )
        r_i = np.tanh(ARGS.score_ingredients_level_mag 
                      * (curPrice.ingredientsLevelForLuxuries - 1))
        if mm > self.ExpectedStandardProfitRate:
            r -= r_c * (1 + r_i)
        else:
            r = r - r_c + \
                np.tanh(ARGS.score_wrong_ingredients_level_mag 
                        * ((curPrice.ingredientsLevelForLuxuries - 1) ** 2))
                     
        r += np.tanh(ARGS.score_profit_rate_mag *
                       ((curPrice.standardProfitRate 
                         - self.ExpectedStandardProfitRate) ** 2))

        self._restoreLcIpp()

        if newCompanyOfNecessaries is not None:
            self.necessariesCompanies.pop()
        if newCompanyOfLuxuries is not None:
            self.luxuriesCompanies.pop()
        if newCompanyOfIngredients is not None:
            self.ingredientsCompanies.pop()

        ## 辞書の作成
        return {
            'score': r,
            'socialDebt': sd,
            'demandOfNecessaries': curDemandOfNecessaries,
            'demandOfLabors': curDemandOfLabors,
            'supplyOfLabors': curSupplyOfLabors,
            'newCompanyOfNecessaries': newCompanyOfNecessaries,
            'productsOfNecessaries': curProductsOfNecessaries,
            'newCompanyOfLuxuries': newCompanyOfLuxuries,
            'productsOfLuxuries': curProductsOfLuxuries,
            'newCompanyOfIngredients': newCompanyOfIngredients,
            'productsOfIngredients': curProductsOfIngredients,
        }


    def step (self, cur_trial='?'):
        minScore = None
        minArg = None

        if ARGS.opt_try == "none":
            res1 = None
        elif ARGS.opt_try == "prev":
            minScore = None
            minArg = None
            for p in self.priceHistory:
                st = self._calcNextState(p.toArray())
                if minScore is None or st['score'] < minScore:
                    minScore = st['score']
                    minArg = p
            res1 = minimize(lambda l: self._calcNextState(l)['score'],
                            minArg.toArray(), method='Nelder-Mead')
        elif ARGS.opt_try == "random":
            minScore = None
            minArg = None
            for i in range(10):
                x = [rand1to100(), rand1to100(), rand1to100(), rand1to100(),
                     self.ExpectedStandardProfitRate, rand1to100()]
                st = self._calcNextState(x)
                if minScore is None or st['score'] < minScore:
                    minScore = st['score']
                    minArg = Price(x)
            res1 = minimize(lambda l: self._calcNextState(l)['score'],
                            minArg.toArray(), method='Nelder-Mead')
        elif ARGS.opt_try == "bh0":
            res1 = scipy.optimize.basinhopping \
                (lambda l: self._calcNextState(l)['score'],
                 self._calcNextInitialValue())
        elif ARGS.opt_try == "de":
            res1 = scipy.optimize.differential_evolution \
                (lambda l: self._calcNextState(l)['score'],
                 [(1, 100)] * 4 + [(0.0, 1.0), (1.0, 100)])
        elif ARGS.opt_try == "shgo":
            res1 = scipy.optimize.shgo \
                (lambda l: self._calcNextState(l)['score'],
                 [(1, 100)] * 4 + [(0.0, 1.0), (1.0, 100)])
        elif ARGS.opt_try == "da":
            res1 = scipy.optimize.dual_annealing \
                (lambda l: self._calcNextState(l)['score'],
                 [(1, 100)] * 4 + [(0.0, 1.0), (1.0, 100)])
        elif ARGS.opt_try == "da0":
            res1 = scipy.optimize.dual_annealing \
                (lambda l: self._calcNextState(l)['score'],
                 [(1, 100)] * 4 + [(0.0, 1.0), (1.0, 100)],
                 x0=self._calcNextInitialValue())
        else:
            raise ValueError("Illegal OPT_TRY: " + ARGS.opt_try)

        res2 = minimize(lambda l: self._calcNextState(l)['score'],
                        self._calcNextInitialValue(clip=False),
                        method='Nelder-Mead')

        if res1 is None:
            res1 = res2

        self.debug = True
        st1 = self._calcNextState(res1.x)
        st2 = self._calcNextState(res2.x)
        if st1['score'] <= st2['score']:
            st = st1
            res = res1
        else:
            st = st2
            res = res2

        if np.isinf(st['score']):
            raise ValueError("Economy Collapsed.")

#        st = self._calcNextState(res.x)
        if not hasattr(res, 'success') or res.success:
            print("[{0}] OptForStep: iterated {1} times score={2}"
                  .format(cur_trial, res.nit, st['score']))
        else:
            print("[{0}] OptForStep(Fail): iterated {1} times score={2}"
                  .format(cur_trial, res.nit, st['score']))
        self.debug = False
        curPrice = Price(res.x)
        socialDebt = 0
        totalSavings = 0
        for x in self.labors:
            totalSavings += x.savings
        for x in self.capitalists:
            totalSavings += x.savings
            totalSavings -= x.debt
        prevTotalSavings = totalSavings
        

        ## 新企業の更新
        if st['newCompanyOfNecessaries'] is not None:
            self.necessariesCompanies.append(st['newCompanyOfNecessaries'])
            self.lastCompanyOfNecessaries = st['newCompanyOfNecessaries']
        if st['newCompanyOfLuxuries'] is not None:
            self.luxuriesCompanies.append(st['newCompanyOfLuxuries'])
            self.lastCompanyOfLuxuries = st['newCompanyOfLuxuries']
        if st['newCompanyOfIngredients'] is not None:
            self.ingredientsCompanies.append(st['newCompanyOfIngredients'])
            self.lastCompanyOfIngredients = st['newCompanyOfIngredients']

        self._setTemporalLcIpp(curPrice.ingredientsLevelForLuxuries,
                               self.workersLevel)

        ## 貯蓄の更新
        trueWages = None
        if st['demandOfLabors'] > st['supplyOfLabors']:
            trueWages = (2 * curPrice.wages \
                             * (st['demandOfLabors'] - st['supplyOfLabors']) \
                             / st['supplyOfLabors']) + curPrice.wages
            socialDebt += 2 * curPrice.wages \
                * (st['demandOfLabors'] - st['supplyOfLabors'])
        else:
            socialDebt += (st['supplyOfLabors'] - st['demandOfLabors']) \
                * curPrice.wages
            trueWages = curPrice.wages
        curD = st['demandOfNecessaries'] \
            / (len(self.capitalists) + len(self.labors))
        curS = curPrice.wages - curPrice.necessaries * curD
        notWorking = len(self.labors) - st['supplyOfLabors']
        base = curS * 3.0
        if base < 0:
            base = 0
        curDemandOfLuxuries = 0 ## デバッグ用
        for i, x in enumerate(self.labors):
            if i < st['supplyOfLabors']:
                x.working = True
                if x.savings + curS >= base \
                        and ((x.savings + curS - base) / 3) + base \
                        >= curPrice.luxuries:
                    curDemandOfLuxuries += 1
                    x.savings += trueWages - curPrice.necessaries * curD \
                        - curPrice.luxuries
                else:
                    x.savings += trueWages - curPrice.necessaries * curD
            else:
                x.working = False
                socialDebt += curPrice.necessaries * curD ## 生活保護
                if (x.savings >= base
                    and ((x.savings - base) / 3) + base >= curPrice.luxuries):
                    curDemandOfLuxuries += 1
                    x.savings -= curPrice.luxuries
        for i, x in enumerate(self.capitalists):
            x.working = True
            # socialDebt += curPrice.wages
            if x.savings + curS >= base \
                    and ((x.savings + curS - base) / 3) + base \
                    >= curPrice.luxuries:
                curDemandOfLuxuries += 1
                x.savings += curPrice.wages - curPrice.necessaries * curD \
                    - curPrice.luxuries
            else:
                x.savings += curPrice.wages - curPrice.necessaries * curD

        if abs(curDemandOfLuxuries - sum(st['productsOfLuxuries'])) > 0.1:
            print("Error about Luxuries!")

        ## 利潤履歴の更新
        print("Companies N:{0}, L:{1}, I:{2}"
              .format(len(self.necessariesCompanies),
                      len(self.luxuriesCompanies),
                      len(self.ingredientsCompanies)))
        print("[Type] [OwnerId] [ProfitRate]")
        curP = st['productsOfNecessaries']
        p = curPrice.necessaries
        l = self.profitHistoryOfNecessaries
        prs = []
        maxPr = None
        for i,c in enumerate(self.necessariesCompanies):
            pr = self._calcProfitOfCompany(c, p, curPrice, curP[i])
            prs.append(pr)
            print("N {0} {1}".format(c.ownerId, pr))
            if maxPr is None or maxPr < pr:
                maxPr = pr
        l.append(maxPr)
        if len(l) > 10:
            l = l[-10:]
        self.profitHistoryOfNecessaries = l
        profitsOfNecessaries = prs

        curP = st['productsOfLuxuries']
        p = curPrice.luxuries
        l = self.profitHistoryOfLuxuries
        prs = []
        maxPr = None
        for i,c in enumerate(self.luxuriesCompanies):
            pr = self._calcProfitOfCompany(c, p, curPrice, curP[i])
            prs.append(pr)
            print("L {0} {1}".format(c.ownerId, pr))
            if maxPr is None or maxPr < pr:
                maxPr = pr
        l.append(maxPr)
        if len(l) > 10:
            l = l[-10:]
        self.profitHistoryOfLuxuries = l
        profitsOfLuxuries = prs

        curP = st['productsOfIngredients']
        p = curPrice.ingredients
        l = self.profitHistoryOfIngredients
        prs = []
        maxPr = None
        for i,c in enumerate(self.ingredientsCompanies):
            pr = self._calcProfitOfCompany(c, p, curPrice, curP[i])
            prs.append(pr)
            print("I {0} {1}".format(c.ownerId, pr))
            if maxPr is None or maxPr < pr:
                maxPr = pr
        l.append(maxPr)
        if len(l) > 10:
            l = l[-5:]
        self.profitHistoryOfIngredients = l
        profitsOfIngredients = prs

        ## 固定資本新設
        newC = st['newCompanyOfNecessaries']
        if (newC is not None):
            self.capitalists[newC.ownerId].debt += newC.fixedAssetsPrice
        newC = st['newCompanyOfLuxuries']
        if (newC is not None):
            self.capitalists[newC.ownerId].debt += newC.fixedAssetsPrice
        newC = st['newCompanyOfIngredients']
        if (newC is not None):
            self.capitalists[newC.ownerId].debt += newC.fixedAssetsPrice

        ## 資本家利潤
        profit = [0 for i in self.capitalists]
        debtable = [0 for i in self.capitalists]
        curP = st['productsOfNecessaries']
        p = curPrice.necessaries
        for i, c in enumerate(self.necessariesCompanies):
            debtable[c.ownerId] += c.fixedAssetsPrice / c.fixedAssetsLife
            profit[c.ownerId] = curP[i] * (\
                p - curPrice.wages * c.laborsPerProduct \
                    + curPrice.ingredients * c.ingredientsPerProduct)

        curP = st['productsOfLuxuries']
        p = curPrice.luxuries
        for i, c in enumerate(self.luxuriesCompanies):
            debtable[c.ownerId] += c.fixedAssetsPrice / c.fixedAssetsLife
            profit[c.ownerId] = curP[i] * (\
                p - curPrice.wages * c.laborsPerProduct \
                    + curPrice.ingredients * c.ingredientsPerProduct)

        curP = st['productsOfIngredients']
        p = curPrice.ingredients
        for i, c in enumerate(self.ingredientsCompanies):
            debtable[c.ownerId] += c.fixedAssetsPrice / c.fixedAssetsLife
            profit[c.ownerId] = curP[i] * (\
                p - curPrice.wages * c.laborsPerProduct \
                    + curPrice.ingredients * c.ingredientsPerProduct)

        for i, cp in enumerate(self.capitalists):
            profit[i] -= curPrice.wages
            debtable[i] += curPrice.wages
            socialDebt += debtable[i]
            cp.debt -= debtable[i]
            if cp.debt < 0:
                socialDebt += cp.debt
                cp.debt = 0
            if profit[i] < 0:
                socialDebt += - profit[i]
            else:
                cp.debt -= profit[i]
                if cp.debt < 0:
                    socialDebt += cp.debt
                    cp.debt = 0

        self._restoreLcIpp()

        ## 固定資本期間終了
        companies = self.necessariesCompanies
        newC = st['newCompanyOfNecessaries']
        l = []
        for c in companies:
            c.fixedAssetsLifeRest -= 1
            if c.fixedAssetsLifeRest >= 0:
                c.superiority *= c.superiorityDecay
                l.append(c)
        self.necessariesCompanies = l

        companies = self.luxuriesCompanies
        newC = st['newCompanyOfLuxuries']
        l = []
        for c in companies:
            c.fixedAssetsLifeRest -= 1
            if c.fixedAssetsLifeRest >= 0:
                c.superiority *= c.superiorityDecay
                l.append(c)
        self.luxuriesCompanies = l

        companies = self.ingredientsCompanies
        newC = st['newCompanyOfIngredients']
        l = []
        for c in companies:
            c.fixedAssetsLifeRest -= 1
            if c.fixedAssetsLifeRest >= 0:
                c.superiority *= c.superiorityDecay
                l.append(c)
        self.ingredientsCompanies = l

        totalSavings = 0
        for x in self.labors:
            totalSavings += x.savings
        for x in self.capitalists:
            totalSavings += x.savings
            totalSavings -= x.debt

        ## レポート
        print("")
        if ARGS.opt_try != "none":
            print("Score {0} vs {1}"
                  .format(st1['score'], st2['score']))
        print("Increase of Social Debt: {0}"
              .format(socialDebt))
        print("Increase of Total Savings: {0}"
              .format(totalSavings - prevTotalSavings))
        print("Price W:{0}, N:{1}, L:{2}, I:{3}, S:{4}, IL:{5}"
              .format(curPrice.wages, curPrice.necessaries,
                      curPrice.luxuries, curPrice.ingredients,
                      curPrice.standardProfitRate,
                      curPrice.ingredientsLevelForLuxuries))
        print("Worker's Savings Increase: {0}".format(curS))
        print("Supply of Labors = {0}".format(st['supplyOfLabors']))
        print("(Demand - Supply) of Labors = {0}"
              .format(st['demandOfLabors'] - st['supplyOfLabors']))
        print("Profit Rate N:{0}, L:{1}, I:{2}"
              .format(self.profitHistoryOfNecessaries[-1],
                      self.profitHistoryOfLuxuries[-1],
                      self.profitHistoryOfIngredients[-1]))
        print("Product N:{0}, L:{1}, I:{2}"
              .format(sum(st['productsOfNecessaries']),
                      sum(st['productsOfLuxuries']),
                      sum(st['productsOfIngredients'])))
        savings = list(map(lambda x: float(x.savings), self.labors \
                               + self.capitalists))
        print("Savings mean:{0}, variance:{1}".format(np.mean(savings),
                                                      np.var(savings)))
        print("", flush=True)


        ## socialDebt の計算間違いがないかチェック。
        if not np.isclose(st['socialDebt'], socialDebt):
            raise ValueError("Social Debt doesn't fit.: {0}:{1}"
                             .format(st['socialDebt'], socialDebt))


        ## self.history に記録。
        h = self.history
        h['score'].append(st['score'])
        if ARGS.opt_try != "none":
            h['score_1'].append(st1['score'])
            h['score_2'].append(st2['score'])

        for l, d in zip(['price_W', 'price_N', 'price_L', 'price_I',
                         'standard_profit_rate',
                         'ingredients_level_for_luxuries'],
                        [curPrice.wages, curPrice.necessaries,
                         curPrice.luxuries, curPrice.ingredients,
                         curPrice.standardProfitRate,
                         curPrice.ingredientsLevelForLuxuries]):
            h[l].append(d)

        for l, d in zip(['num_companies_N', 'num_companies_L',
                         'num_companies_I'],
                        [len(self.necessariesCompanies),
                         len(self.luxuriesCompanies),
                         len(self.ingredientsCompanies)]):
            h[l].append(d)

        h['social_debt'].append(socialDebt)
        h['savings'].append(totalSavings)
        h['savings_increase'].append(totalSavings - prevTotalSavings)
        h['workers_savings_increase'].append(curS)
        h['supply_of_labors'].append(st['supplyOfLabors'])
        h['demand_vs_supply_of_labors'].append(st['demandOfLabors']
                                               - st['supplyOfLabors'])
    
        for l, d in zip(['max_profit_rate_N', 'max_profit_rate_L',
                         'max_profit_rate_I', 'max_profit_rate_mean',
                         'max_profit_rate_all'],
                        [self.profitHistoryOfNecessaries[-1],
                         self.profitHistoryOfLuxuries[-1],
                         self.profitHistoryOfIngredients[-1],
                         np.mean([self.profitHistoryOfNecessaries[-1],
                                  self.profitHistoryOfLuxuries[-1],
                                  self.profitHistoryOfIngredients[-1]]),
                         max(profitsOfNecessaries + profitsOfLuxuries
                             + profitsOfIngredients)]):
            h[l].append(d)

        for l, d in zip(['mean_profit_rate_N', 'mean_profit_rate_L',
                         'mean_profit_rate_I', 'mean_profit_rate_all'],
                        [np.mean(profitsOfNecessaries),
                         np.mean(profitsOfLuxuries),
                         np.mean(profitsOfIngredients),
                         np.mean(profitsOfNecessaries + profitsOfLuxuries
                                 + profitsOfIngredients)]):
            h[l].append(d)

        for l, d in zip(['product_N', 'product_L', 'product_I'],
                        [sum(st['productsOfNecessaries']),
                         sum(st['productsOfLuxuries']),
                         sum(st['productsOfIngredients'])]):
            h[l].append(d)

        h['savings_mean'].append(np.mean(savings))
        h['savings_variance'].append(np.var(savings))
        
        h['workers_level'].append(self.workersLevel)


        ## 労働者と失業者をランダムにソート
        working = []
        notWorking = []
        for x in randomList(self.labors):
            if x.working:
                working.append(x)
            else:
                notWorking.append(x)
        self.labors = working + notWorking

        ## workersLevel の更新
        if all(map(lambda x: len(x) >= ARGS.expected_companies,
                   [self.necessariesCompanies,
                    self.luxuriesCompanies,
                    self.ingredientsCompanies])):
            self.workersLevel *= ARGS.update_workers_level

        ## 履歴の更新
        self.prevPrice = curPrice
        self.priceHistory.append(curPrice)
        self.prevDemandOfNecessaries = st['demandOfNecessaries']
        self.prevSocialDebt = socialDebt
        self.cumulativeSocialDebt += socialDebt

        ## 乱数の更新
        self.randomstate = random.getstate()


    def initialize_0 (self):
        ## 0番目に決めた初期値の条件:
        ## 労働者の人口は 1000。資本家は 5人で別々の戦略。
        ## 労働者の賃金 30、生活費 10、貯蓄は平均 10 の指数分布。
        ## 前期の必需品需要は 2010
        ## 前期の価格は賃金 30、必需品 10、贅沢品 15、原料 15 とする。
        ## 失業者は 1000 人のうち 300人。
        ## 会社は各1社で、オーナーはランダムで、LifeRest もランダム。
        ## 利潤率はどの会社も前期・前々期と「市場利子率」と同じだったとする。
        ## 固定資産の原料は 5 単位で価格は 5 * 5 = 25。
        ## 優位性減衰は 0.9 で、優位性は 0.3 を期間分減衰したもの。
        ## 商品1個あたりの原料は 0.5、商品1個あたりの労働は 0.3 から、
        ## 目やすに「いい感じになるよう少しずつ動かした。

        self.cumulativeSocialDebt = 0
        self.prevSocialDebt = 0

        ## 労働者は1000人。失業者は300人。
        numLabors = 1000
        numNotWorking = 300
        ## 前期の価格 p と必需品需要を適当に…。
        p = Price([30, 10, 15, 15, self.ExpectedStandardProfitRate, 1.0])
        self.prevPrice = p
        self.priceHistory.append(p)
        self.prevDemandOfNecessaries = 2010
        ## 貯蓄 s は前期新規貯蓄。
        s = p.wages - p.necessaries * self.prevDemandOfNecessaries \
            / (numLabors + 5)

        ## 資本家は 5 人。
        self.capitalists = []

        cp = Capitalist()
        self.capitalists.append(cp)
        cp.id = 0
        cp.savings = random.expovariate(1.0/s)
        cp.debt = 0
        cp.makeCompany = cp._makeCompany_r

        cp = Capitalist()
        self.capitalists.append(cp)
        cp.id = 1
        cp.savings = random.expovariate(1.0/s)
        cp.debt = 0
        cp.makeCompany = cp._makeCompany_0_0

        cp = Capitalist()
        self.capitalists.append(cp)
        cp.id = 2
        cp.savings = random.expovariate(1.0/s)
        cp.debt = 0
        cp.makeCompany = cp._makeCompany_0_1

        cp = Capitalist()
        self.capitalists.append(cp)
        cp.id = 3
        cp.savings = random.expovariate(1.0/s)
        cp.debt = 0
        cp.makeCompany = cp._makeCompany_1_0

        cp = Capitalist()
        self.capitalists.append(cp)
        cp.id = 4
        cp.savings = random.expovariate(1.0/s)
        cp.debt = 0
        cp.makeCompany = cp._makeCompany_1_1

        ## 労働者を1000人創り、貯蓄を決定。
        self.labors = []
        for i in range(numLabors):
            x = Labor()
            self.labors.append(x)
            if i >= numLabors - numNotWorking:
                x.working = False
            else:
                x.worknig = True
            x.savings = random.expovariate(1.0/s)

        ## 利潤率はどの会社も「市場利子率」だったとする。
        r = self.ExpectedStandardProfitRate
        self.profitHistoryOfNecessaries = [r, r]
        self.profitHistoryOfLuxuries = [r, r]
        self.profitHistoryOfIngredients = [r, r]
            
        ## 会社の条件は同じ。
        self.necessariesCompanies = []
        c = Company()
        self.necessariesCompanies.append(c)
        self.lastCompanyOfNecessaries = c
        c.ownerId = random.randint(0, 4)
        c.fixedAssetsIngredients = 5
        c.fixedAssetsPrice = c.fixedAssetsIngredients * p.ingredients
        c.fixedAssetsLife = 10
        c.fixedAssetsLifeRest = random.randint(0, 9)
        c.superiorityDecay = 0.9
        c.superiority = 0.3 * (0.9 ** (c.fixedAssetsLife \
                                           - c.fixedAssetsLifeRest))
        c.laborsPerProduct = 0.1 * ARGS.init_lpp_mag
        c.ingredientsPerProduct = 0.3

        self.luxuriesCompanies = []
        c = Company()
        self.luxuriesCompanies.append(c)
        self.lastCompanyOfLuxuries = c
        c.ownerId = random.randint(0, 4)
        c.fixedAssetsIngredients = 5
        c.fixedAssetsPrice = c.fixedAssetsIngredients * p.ingredients
        c.fixedAssetsLife = 10
        c.fixedAssetsLifeRest = random.randint(0, 9)
        c.superiorityDecay = 0.9
        c.superiority = 0.3 * (0.9 ** (c.fixedAssetsLife \
                                           - c.fixedAssetsLifeRest))
        c.laborsPerProduct = 0.3 * ARGS.init_lpp_mag
        c.ingredientsPerProduct = 0.5

        self.ingredientsCompanies = []
        c = Company()
        self.ingredientsCompanies.append(c)
        self.lastCompanyOfIngredients = c
        c.ownerId = random.randint(0, 4)
        c.fixedAssetsIngredients = 5
        c.fixedAssetsPrice = c.fixedAssetsIngredients * p.ingredients
        c.fixedAssetsLife = 10
        c.fixedAssetsLifeRest = random.randint(0, 9)
        c.superiorityDecay = 0.9
        c.superiority = 0.3 * (0.9 ** (c.fixedAssetsLife \
                                           - c.fixedAssetsLifeRest))
        c.laborsPerProduct = 0.3 * ARGS.init_lpp_mag
        c.ingredientsPerProduct = 0.3

        self.randomstate = random.getstate()


if __name__ == '__main__':
    economy = Economy()
    economy.initialize_0()
    print("Acronyms: W: Wages, N: Necessaries, L: Luxuries, I: Ingredients, S: Standard Profit Rate.\n")
    try:
        for i in range(ARGS.trials):
            economy.step(cur_trial=(i+1))
    except:
        traceback.print_exc()

    if ARGS.save_history is not None:
        save_history(ARGS.save_history, economy.history)


