#!/usr/bin/python3
__version__ = '0.0.7' # Time-stamp: <2020-02-12T12:25:42Z>
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

import math
import random
from scipy.optimize import minimize, minimize_scalar
from statistics import mean, variance

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
        if igLimit != None and r > igLimit:
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

    def toArray (self):
        return [self.wages, self.necessaries, self.luxuries, self.ingredients]

def randomList (l):
    return random.sample(l, len(l))

class Economy:
    necessariesElasticity = - 0.6
    surplusVsNecessariesRatio = 0.5
    standardProfitRate = 0.06
    wagesElasticPower = 1.0
    ingredientsLimit = 0.9

    def __init__ (self):
        random.seed()
        self.debug = False
        self.randomstate = random.getstate()
        self.socialDebt = 0
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

    def _calcProfitOfCompany (self, company, price, curPrice, num):
        c = company
        cost = (c.fixedAssetsPrice / c.fixedAssetsLife) \
            + (curPrice.wages * c.laborsPerProduct
               + curPrice.ingredients * c.ingredientsPerProduct) * num
        return (price * num - cost) / cost

    def _calcProducts (self, curPrice, price, demand, companies):
        curS = [0 for c in companies]
        delta = Economy.wagesElasticPower \
            * (curPrice.wages - self.prevPrice.wages) / self.prevPrice.wages
        for i, c in enumerate(companies):
            if (price - curPrice.wages * c.laborsPerProduct \
                    - curPrice.ingredients * c.ingredientsPerProduct) > 0:
                curS[i] = c.superiority * (c.laborsPerProduct ** (- delta))
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

    def _calcNextState (self, priceArray):
        random.setstate(self.randomstate)

        curPrice = Price(priceArray)

        if any(map(lambda x: x < 1 or x > 1000, priceArray)):
            return {'score': float('inf')}

        ## 必需品需要の計算
        prevD = self.prevDemandOfNecessaries \
            / (len(self.capitalists) + len(self.labors))
        curD = (1 + ((curPrice.necessaries - self.prevPrice.necessaries)
                     / self.prevPrice.necessaries)
                * Economy.necessariesElasticity) * prevD
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
            q = ((curS - prevS) / prevS) * Economy.surplusVsNecessariesRatio \
                + ((curPrice.necessaries - self.prevPrice.necessaries) \
                       / self.prevPrice.necessaries) \
                       * (1 - Economy.surplusVsNecessariesRatio)
            r = 0.5 * 2 * math.atan(2 * q) / math.pi
            working = 0
            for l in self.labors:
                if l.working:
                    working += 1
            notWorking = len(self.labors) - working
            if notWorking > working:
                curSupplyOfLabors = working + math.ceil(r * working)
            else:
                curSupplyOfLabors = len(self.labors) \
                    - (notWorking - math.floor(r * notWorking))

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
            newC = cp.makeCompany(lastC, curPrice, None)
            newA = companies + [newC]
            curP = self._calcProducts(curPrice, p, u, newA)
            ## その期の利潤率と3期平均利潤率の低い方が市場利子率より
            ## 高いとき投資する。
            pr = self._calcProfitOfCompany(newC, p, curPrice, curP[-1])
            pr2 = sum(prHistory[-2:] + [pr]) / 3
            if pr > pr2:
                pr = pr2
            ## 会社が 0 の場合は最後の資本家が責任を持つ。
            if  pr >= Economy.standardProfitRate \
                    or (len(companies) == 0 \
                            and i == len(self.capitalists) - 1):
                companies.append(newC)
                break
            else:
                newC = None
        if newC == None:
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
            newC = cp.makeCompany(lastC, curPrice, None)
            newA = companies + [newC]
            curP = self._calcProducts(curPrice, p, u, newA)
            ## その期の利潤率と3期平均利潤率の低い方が市場利子率より
            ## 高いとき投資する。
            pr = self._calcProfitOfCompany(newC, p, curPrice, curP[-1])
            pr2 = sum(prHistory[-2:] + [pr]) / 3
            if pr > pr2:
                pr = pr2
            ## 会社が 0 の場合は最後の資本家が責任を持つ。
            if  pr >= Economy.standardProfitRate \
                    or (len(companies) == 0 \
                            and i == len(self.capitalists) - 1):
                companies.append(newC)
                break
            else:
                newC = None
        if newC == None:
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
        if newCompanyOfNecessaries != None:
            u += newCompanyOfNecessaries.fixedAssetsIngredients
        if newCompanyOfLuxuries != None:
            u += newCompanyOfLuxuries.fixedAssetsIngredients
        for i, cp in enumerate(sorted(randomList(self.capitalists),
                                      key=self._cp_debt_key)):
            ## 投資のためには debt がないのが条件
            if cp.debt >= 0.1 \
                    and not (len(companies) == 0 \
                                 and i == len(self.capitalists) - 1):
                continue
            newC = cp.makeCompany(lastC, curPrice, Economy.ingredientsLimit)
            newA = companies + [newC]
            curP = self._calcProductsOfIngredients(\
                curPrice, p, u + newC.fixedAssetsIngredients, newA)
            ## その期の利潤率と3期平均利潤率の低い方が市場利子率より
            ## 高いとき投資する。
            pr = self._calcProfitOfCompany(newC, p, curPrice, curP[-1])
            pr2 = sum(prHistory[-2:] + [pr]) / 3
            if pr > pr2:
                pr = pr2
            ## 会社が 0 の場合は最後の資本家が責任を持つ。
            if  pr >= Economy.standardProfitRate \
                    or (len(companies) == 0 \
                            and i == len(self.capitalists) - 1):
                companies.append(newC)
                break
            else:
                newC = None
        if newC == None:
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

        ## スコアの計算
        r = 0
        ## チューンアップ前:
        r += (curDemandOfLabors - curSupplyOfLabors) ** 2
        ## チューンアップ後:
        ## 3 倍にする意味はあまり明らかではない。
#        r += 3 * ((curDemandOfLabors - curSupplyOfLabors) ** 2)

        ## スコアを socialDebt に近いものにする実験をしたことがあった。
#        if curDemandOfLabors > curSupplyOfLabors:
#            r += 2 * curPrice.wages * (curDemandOfLabors - curSupplyOfLabors)
#        else:
#            r += (curSupplyOfLabors - curDemandOfLabors) * curPrice.wages

        ## 新規貯蓄をスコアに加える実験をしたこともあった。
#        r -= curSavingsOfLabors * curSupplyOfLabors

        ## 最適化するのに利潤総額を見るべきか、利益率を見るべきか…。
        ## 本来は利益率は起業・退出によって平準化されるべきである。しか
        ## し、そうすると利益率がマイナスのまま、参入のおきない分野が出
        ## てくる…。なぜだ？
        pr = None
        curP = curProductsOfNecessaries
        p = curPrice.necessaries
        r1 = 0
        for i, c in enumerate(self.necessariesCompanies):
            pr = - (c.fixedAssetsPrice / c.fixedAssetsLife) \
                + (p - c.laborsPerProduct * curPrice.wages \
                       -  c.ingredientsPerProduct * curPrice.ingredients) \
                       * curP[i]
            r1 += pr
        ## チューンアップ前:
        r -= r1
        ## チューンアップ後:
#        if r1 < 0:
#            r -= 10 * r1 / len(self.necessariesCompanies)
#        else:
#            r -= r1 / len(self.necessariesCompanies)

        curP = curProductsOfLuxuries
        p = curPrice.luxuries
        r1 = 0
        for i, c in enumerate(self.luxuriesCompanies):
            pr = - (c.fixedAssetsPrice / c.fixedAssetsLife) \
                + (p - c.laborsPerProduct * curPrice.wages \
                       -  c.ingredientsPerProduct * curPrice.ingredients) \
                       * curP[i]
            r1 += pr
        ## チューンアップ前:
        r -= r1
        ## チューンアップ後:
#        if r1 < 0:
#            r -= 10 * r1 / len(self.luxuriesCompanies)
#        else:
#            r -= r1 / len(self.luxuriesCompanies)

        curP = curProductsOfIngredients
        p = curPrice.ingredients
        r1 = 0
        for i, c in enumerate(self.ingredientsCompanies):
            pr = - (c.fixedAssetsPrice / c.fixedAssetsLife) \
                + (p - c.laborsPerProduct * curPrice.wages \
                       -  c.ingredientsPerProduct * curPrice.ingredients) \
                       * curP[i]
            r1 += pr
        ## チューンアップ前:
        r -= r1
        ## チューンアップ後:
#        if r1 < 0:
#            r -= 10 * r1 / len(self.ingredientsCompanies)
#        else:
#            r -= r1 / len(self.ingredientsCompanies)

        if newCompanyOfNecessaries != None:
            self.necessariesCompanies.pop()
        if newCompanyOfLuxuries != None:
            self.luxuriesCompanies.pop()
        if newCompanyOfIngredients != None:
            self.ingredientsCompanies.pop()

        ## 辞書の作成
        return {
            'score': r,
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


    def step (self):
        res = minimize(lambda l: self._calcNextState(l)['score'],
                       self.prevPrice.toArray(), method='Nelder-Mead')
        self.debug = True
        st = self._calcNextState(res.x)
        if res.success:
            print("OptForStep: iterated {0} times score={1}"
                  .format(res.nit, st['score']))
        else:
            print("OptForStep(Fail): iterated {0} times score={1}"
                  .format(res.nit, st['score']))
        self.debug = False
        curPrice = Price(res.x)
        prevSocialDebt = self.socialDebt
        totalSavings = 0
        for x in self.labors:
            totalSavings += x.savings
        for x in self.capitalists:
            totalSavings += x.savings
            totalSavings -= x.debt
        prevTotalSavings = totalSavings


        ## 新企業の更新
        if st['newCompanyOfNecessaries'] != None:
            self.necessariesCompanies.append(st['newCompanyOfNecessaries'])
            self.lastCompanyOfNecessaries = st['newCompanyOfNecessaries']
        if st['newCompanyOfLuxuries'] != None:
            self.luxuriesCompanies.append(st['newCompanyOfLuxuries'])
            self.lastCompanyOfLuxuries = st['newCompanyOfLuxuries']
        if st['newCompanyOfIngredients'] != None:
            self.ingredientsCompanies.append(st['newCompanyOfIngredients'])
            self.lastCompanyOfIngredients = st['newCompanyOfIngredients']

        ## 貯蓄の更新
        trueWages = None
        if st['demandOfLabors'] > st['supplyOfLabors']:
            trueWages = (2 * curPrice.wages \
                             * (st['demandOfLabors'] - st['supplyOfLabors']) \
                             / st['supplyOfLabors']) + curPrice.wages
            self.socialDebt += 2 * curPrice.wages \
                * (st['demandOfLabors'] - st['supplyOfLabors'])
        else:
            self.socialDebt += (st['supplyOfLabors'] - st['demandOfLabors']) \
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
                self.socialDebt += curPrice.necessaries * curD ## 生活保護
                if (x.savings >= base
                    and ((x.savings - base) / 3) + base >= curPrice.luxuries):
                    curDemandOfLuxuries += 1
                    x.savings -= curPrice.luxuries
        for i, x in enumerate(self.capitalists):
            x.working = True
            # self.socialDebt += curPrice.wages
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
        maxPr = None
        for i,c in enumerate(self.necessariesCompanies):
            pr = self._calcProfitOfCompany(c, p, curPrice, curP[i])
            print("N {0} {1}".format(c.ownerId, pr))
            if maxPr == None or maxPr < pr:
                maxPr = pr
        l.append(maxPr)
        if len(l) > 10:
            l = l[-10:]
        self.profitHistoryOfNecessaries = l

        curP = st['productsOfLuxuries']
        p = curPrice.luxuries
        l = self.profitHistoryOfLuxuries
        maxPr = None
        for i,c in enumerate(self.luxuriesCompanies):
            pr = self._calcProfitOfCompany(c, p, curPrice, curP[i])
            print("L {0} {1}".format(c.ownerId, pr))
            if maxPr == None or maxPr < pr:
                maxPr = pr
        l.append(maxPr)
        if len(l) > 10:
            l = l[-10:]
        self.profitHistoryOfLuxuries = l

        curP = st['productsOfIngredients']
        p = curPrice.ingredients
        l = self.profitHistoryOfIngredients
        maxPr = None
        for i,c in enumerate(self.ingredientsCompanies):
            pr = self._calcProfitOfCompany(c, p, curPrice, curP[i])
            print("I {0} {1}".format(c.ownerId, pr))
            if maxPr == None or maxPr < pr:
                maxPr = pr
        l.append(maxPr)
        if len(l) > 10:
            l = l[-5:]
        self.profitHistoryOfIngredients = l

        ## 固定資本新設
        newC = st['newCompanyOfNecessaries']
        if (newC != None):
            self.capitalists[newC.ownerId].debt += newC.fixedAssetsPrice
        newC = st['newCompanyOfLuxuries']
        if (newC != None):
            self.capitalists[newC.ownerId].debt += newC.fixedAssetsPrice
        newC = st['newCompanyOfIngredients']
        if (newC != None):
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
            self.socialDebt += debtable[i]
            cp.debt -= debtable[i]
            if cp.debt < 0:
                self.socialDebt += cp.debt
                cp.debt = 0
            if profit[i] < 0:
                self.socialDebt += - profit[i]
            else:
                cp.debt -= profit[i]
                if cp.debt < 0:
                    self.socialDebt += cp.debt
                    cp.debt = 0

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
        print("Increase of Social Debt: {0}"
              .format(self.socialDebt - prevSocialDebt))
        print("Increase of Total Savings: {0}"
              .format(totalSavings - prevTotalSavings))
        print("Price W:{0}, N:{1}, L:{2}, I:{3}"
              .format(curPrice.wages, curPrice.necessaries,
                      curPrice.luxuries, curPrice.ingredients))
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
        savings = map(lambda x: float(x.savings), self.labors \
                          + self.capitalists)
        r = [c for c in savings]  ## なぜこれが必要なの？
        print("Savings mean:{0}, variance:{1}".format(mean(r), variance(r)))
        print("", flush=True)


        ## 労働者と失業者をランダムにソート
        working = []
        notWorking = []
        for x in randomList(self.labors):
            if x.working:
                working.append(x)
            else:
                notWorking.append(x)
        self.labors = working + notWorking

        ## 履歴の更新
        self.prevPrice = curPrice
        self.prevDemandOfNecessaries = st['demandOfNecessaries']

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

        self.socialDebt = 0

        ## 労働者は1000人。失業者は300人。
        numLabors = 1000
        numNotWorking = 300
        ## 前期の価格 p と必需品需要を適当に…。
        p = Price([30, 10, 15, 15])
        self.prevPrice = p
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
        r = Economy.standardProfitRate
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
        c.laborsPerProduct = 0.1
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
        c.laborsPerProduct = 0.3
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
        c.laborsPerProduct = 0.3
        c.ingredientsPerProduct = 0.3

        self.randomstate = random.getstate()


if __name__ == '__main__':
    economy = Economy()
    economy.initialize_0()
    print("Acronyms: W: Wages, N: Necessaries, L: Luxuries, I: Ingredients.\n")
    for i in range(20):
        economy.step()
