{-
 ----------------------------------------------------------------------------------
 -  Copyright (C) 2010-2011  Massachusetts Institute of Technology
 -  Copyright (C) 2010-2011  Yuan Tang <yuantang@csail.mit.edu>
 -                           Charles E. Leiserson <cel@mit.edu>
 -   
 -   This program is free software: you can redistribute it and/or modify
 -   it under the terms of the GNU General Public License as published by
 -   the Free Software Foundation, either version 3 of the License, or
 -   (at your option) any later version.
 -
 -   This program is distributed in the hope that it will be useful,
 -   but WITHOUT ANY WARRANTY; without even the implied warranty of
 -   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 -   GNU General Public License for more details.
 -
 -   You should have received a copy of the GNU General Public License
 -   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 -
 -   Suggestsions:                  yuantang@csail.mit.edu
 -   Bugs:                          yuantang@csail.mit.edu
 -
 --------------------------------------------------------------------------------
 -}

module PGenUtils where

import Text.ParserCombinators.Parsec
import Control.Monad
import Data.List
import Data.Bits
import PGenData
import qualified Data.Map as Map

updatePTile :: PTile -> ParserState -> ParserState
updatePTile l_tile parserState =
    parserState { pTile = Map.insert (tName l_tile) l_tile $ pTile parserState }

updatePTileOrder :: Int -> ParserState -> ParserState
updatePTileOrder l_tile_order parserState =
    parserState { pTileOrder = l_tile_order }

updatePGuardOrder :: Int -> ParserState -> ParserState
updatePGuardOrder l_guardOrder parserState =
    parserState { pGuardOrder = l_guardOrder }

updatePGuardFuncOrder :: Int -> ParserState -> ParserState
updatePGuardFuncOrder l_guardFuncOrder parserState =
    parserState { pGuardFuncOrder = l_guardFuncOrder }

updatePGenPlanOrder :: Int -> ParserState -> ParserState
updatePGenPlanOrder l_gen_plan_order parserState =
    parserState { pGenPlanOrder = l_gen_plan_order }

updatePKernelFunc :: PKernelFunc -> ParserState -> ParserState
updatePKernelFunc l_kernelFunc parserState =
    parserState { pKernelFunc = Map.insert (kfName l_kernelFunc) l_kernelFunc $ pKernelFunc parserState }

updatePKernel :: PKernel -> ParserState -> ParserState
updatePKernel l_kernel parserState =
    parserState { pKernel = Map.insert (kName l_kernel) l_kernel $ pKernel parserState }

updatePGuardFunc :: PGuardFunc -> ParserState -> ParserState
updatePGuardFunc l_gFunc parserState =
    parserState { pGuardFunc = Map.insert (gfName l_gFunc) l_gFunc $ pGuardFunc parserState }

updatePGuard :: PGuard -> ParserState -> ParserState
updatePGuard l_guard parserState =
    parserState { pGuard = Map.insert (gName l_guard) l_guard $ pGuard parserState }

updateObase :: PMode -> ParserState -> ParserState
updateObase mode parserState = parserState { pMode = mode }

updatePShape :: PShape -> ParserState -> ParserState
updatePShape l_pShape parserState = 
    parserState { pShape = Map.insert (shapeName l_pShape) l_pShape (pShape parserState) }

updatePArray :: [(PName, PArray)] -> ParserState -> ParserState
updatePArray [] parserState = parserState
updatePArray pL@(p:ps) parserState =
    parserState { pArray = foldr pMapInsert (pArray parserState) pL }

updatePStencil :: [(PName, PStencil)] -> ParserState -> ParserState
updatePStencil [] parserState = parserState
updatePStencil pL@(p:ps) parserState =
    parserState { pStencil = foldr pMapInsert (pStencil parserState) pL }

updatePGenPlan :: (Int, PStencil) -> ParserState -> ParserState
updatePGenPlan p parserState =
    parserState { pGenPlan = pMapInsert p (pGenPlan parserState) }

updatePRange :: [(PName, PRange)] -> ParserState -> ParserState
updatePRange [] parserState = parserState
updatePRange pL@(p:ps) parserState =
    parserState { pRange = foldr pMapInsert (pRange parserState) pL }

pMapInsert :: (Ord k) => (k, a) -> Map.Map k a -> Map.Map k a
pMapInsert (l_key, l_value) l_map = Map.insert l_key l_value l_map

updateStencilArray :: String -> [PArray] -> GenParser Char ParserState String
updateStencilArray _ [] = return ("")
updateStencilArray l_id l_pArrays@(a:as) =
    do updateState $ updateStencilArrayItem l_id a
       updateStencilArray l_id as

updateStencilArrayItem :: String -> PArray -> ParserState -> ParserState
updateStencilArrayItem l_id l_pArray parserState =
    let f k x =
            if sName x == l_id
                then Just $ x { sArrayInUse = union [l_pArray] (sArrayInUse x) }
                else Nothing
    in  parserState { pStencil = Map.updateWithKey f l_id $ pStencil parserState }

updateStencilBoundary :: String -> Bool -> ParserState -> ParserState
updateStencilBoundary l_id l_regBound parserState =
    let f k x =
            if sName x == l_id
                then Just $ x { sRegBound = l_regBound }
                else Nothing
    in  parserState { pStencil = Map.updateWithKey f l_id $ pStencil parserState }

updateTileOrigGuard :: String -> PGuard -> ParserState -> ParserState
updateTileOrigGuard l_id l_guard parserState =
    let f k x =
            if tName x == l_id
                then Just $ x { tOrigGuard = l_guard }
                else Nothing
    in  parserState { pTile = Map.updateWithKey f l_id $ pTile parserState }

updateStencilRegTileKernel :: String -> [(PGuard, PTile)] -> ParserState -> ParserState
updateStencilRegTileKernel l_id l_regTileKernel parserState =
    let f k x =
            if sName x == l_id
                then Just $ x { sRegTileKernel = l_regTileKernel }
                else Nothing
    in  parserState { pStencil = Map.updateWithKey f l_id $ pStencil parserState }

updateStencilToggle :: String -> Int -> ParserState -> ParserState
updateStencilToggle l_id l_toggle parserState =
    let f k x =
            if sName x == l_id
                then Just $ x { sToggle = max l_toggle $ sToggle x }
                else Nothing
    in  parserState { pStencil = Map.updateWithKey f l_id $ pStencil parserState }

updateStencilTimeShift :: String -> Int -> ParserState -> ParserState
updateStencilTimeShift l_id l_timeShift parserState =
    let f k x =
            if sName x == l_id
                then Just $ x { sTimeShift = min l_timeShift $ sTimeShift x }
                else Nothing
    in  parserState { pStencil = Map.updateWithKey f l_id $ pStencil parserState }

updateStencilShape :: String -> PShape -> ParserState -> ParserState
updateStencilShape l_id l_pShape parserState =
    let f k x =
            if sName x == l_id
                then Just $ x { sShape = l_pShape }
                else Nothing
    in  parserState { pStencil = Map.updateWithKey f l_id $ pStencil parserState }

updateArrayBoundary :: String -> Bool -> ParserState -> ParserState
updateArrayBoundary l_id l_regBound parserState =
    let f k x =
            if aName x == l_id
                then Just $ x { aRegBound = l_regBound }
                else Nothing
    in parserState { pArray = Map.updateWithKey f l_id $ pArray parserState }

getToggleFromShape :: [[Int]] -> Int
getToggleFromShape [] = 0
getToggleFromShape l_shapes =
    let l_t = map head l_shapes
        l_t_max = maximum l_t
        l_t_min = minimum l_t
    in  (1 + l_t_max - l_t_min)

getSlopesFromShape :: Int -> [[Int]] -> [Int]
getSlopesFromShape l_height l_shapes = 
    let l_spatials = transpose $ map tail l_shapes
        l_xs = getXSFromShape l_spatials
    in  map (flip div l_height) l_xs

getTimeShiftFromShape :: [[Int]] -> Int
getTimeShiftFromShape [] = 0
getTimeShiftFromShape l_shapes =
    let l_t = map head l_shapes
    in  minimum l_t

getXSFromShape :: [[Int]] -> [Int]
getXSFromShape [] = []
getXSFromShape (a:as) = (maximum $ map abs a):getXSFromShape as

getArrayRegBound :: ParserState -> PArray -> Bool
getArrayRegBound l_state l_pArray =
    case Map.lookup (aName l_pArray) $ pArray l_state of
        Nothing -> False
        Just l_array -> aRegBound l_array
    
getValidKernel :: ParserState -> String -> PKernel
getValidKernel l_state l_kernelName =
    case Map.lookup l_kernelName $ pKernel l_state of
        Nothing -> emptyKernel { kName = l_kernelName, kComment = cUnknown "Kernel" }
        Just l_pKernel -> l_pKernel { kComment = cKnown "Kernel" }

getValidGuard :: ParserState -> String -> PGuard
getValidGuard l_state l_guardName =
    case Map.lookup l_guardName $ pGuard l_state of
        Nothing -> emptyGuard { gName = l_guardName, gComment = [cUnknown "Guard"] }
        Just l_pGuard -> l_pGuard { gComment = [cKnown "Guard"] }

getValidTile :: ParserState -> String -> PTile
getValidTile l_state l_tileName =
    case Map.lookup l_tileName $ pTile l_state of
        Nothing -> 
            let l_kernel = case Map.lookup l_tileName $ pKernel l_state of
                                    Nothing -> emptyKernel
                                    Just l_kern -> l_kern
                l_tileKernel = SK l_kernel
            in  PTile { tName = l_tileName, tRank = kRank l_kernel, tSizes = [1], tKernels = [l_kernel], tComment = "", tOp = PNULL, tOrigGuard = emptyGuard, tOrder = 0, tColor = emptyColor }
        Just l_pTile -> l_pTile { tComment = cKnown "Tile" }

getSetKernelsFromTile :: PTile -> [PKernel]
getSetKernelsFromTile l_tile =
    let l_kernels = tKernels l_tile
    in  map (pSetKernel l_tile) l_kernels
            where pSetKernel l_tile l_kernel = 
                    let l_kfunc = kFunc l_kernel
                        l_gfunc = gFunc $ tOrigGuard l_tile
                        l_tile_op = tOp l_tile
                        l_tile_order = tOrder l_tile
                        l_tile_sizes = tSizes l_tile
                        -- l_res_dims = (length . kfParams) l_kfunc - length l_tile_sizes
                        -- l_res = replicate l_res_dims 0
                        l_res = []
                    in  l_kernel { kFunc = l_kfunc { kfGuardFunc = l_gfunc, kfTileOp = l_tile_op, kfTileOrder = l_tile_order, kfTileSizes = l_tile_sizes ++ l_res, kfTileIndex = kIndex l_kernel ++ l_res } }

getTileKernels :: PTile -> PTileKernel -> [PKernel]
getTileKernels l_tile l_tile_kernels = getKernelsTile [] 0 l_tile l_tile_kernels

getKernelsTile :: [Int] -> Int -> PTile -> PTileKernel -> [PKernel]
getKernelsTile l_indices l_index l_tile (SK l_kernel) =
    let ll_indices = l_indices ++ [l_index]
        ll_kernel_func = kFunc l_kernel
        ll_indices' = if null ll_indices then [] else init ll_indices
        l_gfunc = gFunc $ tOrigGuard l_tile
        l_tile_op = tOp l_tile
        l_tile_order = tOrder l_tile
        ll_kernel = l_kernel { kIndex = ll_indices', kFunc = ll_kernel_func { kfGuardFunc = l_gfunc, kfTileOp = l_tile_op, kfTileOrder = l_tile_order } }
    in  [ll_kernel]
getKernelsTile l_indices l_index l_tile (LK l_tKs@(t:ts)) =
    let ll_indices = l_indices ++ [l_index]
    in  getKernelsTile ll_indices 0 l_tile t ++ 
        getKernelsTile l_indices (l_index + 1) l_tile (LK ts)
getKernelsTile l_indices l_index l_tile (LK []) = []

getTileSizes :: [[Int]] -> [Int]
getTileSizes l_tile_indices = map ((+) 1 . maximum) $ transpose l_tile_indices

getTileDimSizes :: [PKernel] -> ([Int], [PKernel])
getTileDimSizes l_kernels =
    let l_tile_indices = map kIndex l_kernels
        l_tile_sizes = getTileSizes l_tile_indices
    in  (l_tile_sizes, map (pFillTileSizesInPKernel l_tile_sizes) l_kernels)

pFillTileSizesInPKernel :: [Int] -> PKernel -> PKernel
pFillTileSizesInPKernel l_tile_sizes l_kernel = l_kernel { kTileSizes = l_tile_sizes } 

eqTileOpPKernelFunc :: PKernelFunc -> PKernelFunc -> Bool
eqTileOpPKernelFunc a b = (kfTileOp a) == (kfTileOp b)

eqTPKernel :: PKernel -> PKernel -> Bool
eqTPKernel a b = (head $ kIndex a) == (head $ kIndex b)

eqTGroupPKernel :: [PKernel] -> [PKernel] -> Bool
eqTGroupPKernel [] [] = True
eqTGroupPKernel _ [] = False
eqTGroupPKernel [] _ = False
eqTGroupPKernel a b = (head . kIndex . head) a == (head . kIndex . head) b

eqTPTile :: [Int] -> [Int] -> Bool
eqTPTile [] [] = True
eqTPTile [] b = False
eqTPTile a [] = False
eqTPTile a b = head a == head b

eqTileIndex :: [Int] -> [Int] -> Bool
eqTileIndex [] [] = True
eqTileIndex _ [] = True
eqTileIndex [] _ = True
eqTileIndex aL@(a:as) bL@(b:bs) = 
    if a == b
       then eqTileIndex as bs
       else False

pEqList :: Eq a => [a] -> [a] -> PPred
pEqList aL bL = 
    let len_a = length aL
        len_b = length bL
    in  if (len_a < len_b)
           then if aL == take len_a bL then PLEQ else PNEQ
           else if len_a > len_b
                   then if take len_b aL == bL then PGEQ else PNEQ
                   else if aL == bL then PEQ else PNEQ

-- Get the longest common prefix of two lists
pLongestCommonPrefix :: Eq a => [a] -> [a] -> [a]
pLongestCommonPrefix [] [] = []
pLongestCommonPrefix [] _ = []
pLongestCommonPrefix _ [] = []
pLongestCommonPrefix aL@(a:as) bL@(b:bs) = 
    if a == b
       then a:(pLongestCommonPrefix as bs)
       else []

eqKIndex :: PKernel -> PKernel -> Bool
eqKIndex a b = (kIndex a) == (kIndex b)

eqKIndexLen :: PKernel -> PKernel -> Bool
eqKIndexLen a b = (length . kIndex) a == (length . kIndex) b

diffTileOrderPKernel :: PKernel -> PKernel -> Int -> Bool
diffTileOrderPKernel a b n = (abs $ (kfTileOrder . kFunc) a - (kfTileOrder . kFunc) b) == n

eqIndexPKernelTOrder :: PKernel -> PKernel -> Bool
eqIndexPKernelTOrder a b = (eqKIndex a b) && (diffTileOrderPKernel a b 1)

eqIndexPKernel :: PKernel -> PKernel -> Bool
eqIndexPKernel a b = (eqKIndexLen a b) && (diffTileOrderPKernel a b 0)

eqVarLenIndexPKernelTOrder :: PKernel -> PKernel -> Bool
eqVarLenIndexPKernelTOrder a b = (eqTileIndex (kIndex a) (kIndex b)) && (diffTileOrderPKernel a b 1)

pSerializePKernel :: [PKernel] -> Int -> Int -> [PKernel] -> [PKernel]
pSerializePKernel [] n counter l = l
pSerializePKernel kL@(k:ks) n counter l =
    let l_kf = kFunc k
        l_tOrder = (kfTileOrder . kFunc) k
        n' = if l_tOrder == n
                    then n
                    else l_tOrder
        counter' = if l_tOrder == n
                then counter
                else counter+1
        k' = if l_tOrder == n
                then k { kFunc = l_kf { kfTileOrder = counter } }
                else k { kFunc = l_kf { kfTileOrder = counter+1 } }
        l' = l ++ [k']
    in  pSerializePKernel ks n' counter' l'

pFillKIndex :: [Int] -> PKernel -> PKernel
pFillKIndex index k = k { kIndex = index }

pGroupPKernelByMerge :: (PKernel -> PKernel -> Bool) -> [PKernel] -> [[PKernel]]
pGroupPKernelByMerge b [] = []
pGroupPKernelByMerge b l = 
    let l_pKernel_by_tIndex = pGroupPKernelBy eqIndexPKernel l
    in  pGroupPKernelByMergeTerm b $ map (map wrapListMonad) l_pKernel_by_tIndex

wrapListMonad :: a -> [a]
wrapListMonad a = [a]

pGroupPKernelByMergeTerm :: (PKernel -> PKernel -> Bool) -> [[[PKernel]]] -> [[PKernel]]
pGroupPKernelByMergeTerm _ [] = []
pGroupPKernelByMergeTerm _ [a] = a
pGroupPKernelByMergeTerm b ll@(x:y:zs) = 
    -- firstly, let's merge only the following simplified case
    if ((length . kIndex . last . last) x <= (length . kIndex . head . head) y 
        && length x <= length y)
        then let x' = pMergeForward b x y
             in  pGroupPKernelByMergeTerm b (x':zs)
        else if ((length . kIndex . last . last) x > (length . kIndex . head . head) y
               && length x > length y)
                then let x' = pMergeBackward b x y []
                     in  pGroupPKernelByMergeTerm b (x':zs)
                else x ++ pGroupPKernelByMergeTerm b (y:zs)

pMergeForward :: (PKernel -> PKernel -> Bool) -> [[PKernel]] -> [[PKernel]] -> [[PKernel]]
pMergeForward b [] y = y
pMergeForward b x@(a:as) y = 
    let could_merge = foldr (||) False $ 
                        map (b (last a)) (map head y)
    in  if could_merge
           then pMergeForward b as $ pMergeForwardTerm b a y
           else a:(pMergeForward b as y)

pMergeForwardTerm :: (PKernel -> PKernel -> Bool) -> [PKernel] -> [[PKernel]] -> [[PKernel]]
pMergeForwardTerm b x [] = []
pMergeForwardTerm b x yL@(y:ys) =
    if b (last x) (head y)
       then let x' = map (pFillKIndex $ (kIndex . head) y) x
            in  (x' ++ y):(pMergeForwardTerm b x ys) 
       else y:(pMergeForwardTerm b x ys)

pMergeBackward :: (PKernel -> PKernel -> Bool) -> [[PKernel]] -> [[PKernel]] -> [[PKernel]] -> [[PKernel]]
pMergeBackward b x [] l = x ++ l
pMergeBackward b x y@(a:as) l =
    let could_merge = foldr (||) False $
                        map (flip b (head a)) (map last x)
    in  if could_merge 
           then pMergeBackward b (pMergeBackwardTerm b x a) as l
           else pMergeBackward b x as $ l ++ [a]

pMergeBackwardTerm :: (PKernel -> PKernel -> Bool) -> [[PKernel]] -> [PKernel] -> [[PKernel]]
pMergeBackwardTerm b [] y = []
pMergeBackwardTerm b xL@(x:xs) y =
    if b (last x) (head y)
       then let y' = map (pFillKIndex $ (kIndex . last) x) y
            in  (x ++ y'):(pMergeBackwardTerm b xs y)
       else x:(pMergeBackwardTerm b xs y)

pGroupPKernelBy :: (PKernel -> PKernel -> Bool) -> [PKernel] -> [[PKernel]]
pGroupPKernelBy b [] = []
pGroupPKernelBy b l = 
    let n = (kfTileOrder . kFunc . head) l
        l' = pSerializePKernel l n 0 []
    in  pGroupPKernelByTerm b l' []
        
pGroupPKernelByTerm :: (PKernel -> PKernel -> Bool) -> [PKernel] -> [[PKernel]] -> [[PKernel]]
pGroupPKernelByTerm b [] ass = ass
pGroupPKernelByTerm b l@(x:xs) ass =
    let (inSet, outSet) = pPartition b x xs ([x], [])
--  if we use the common partition defined in prelude, because it compares 
--  the tile_order of each PKernel, it will step into an infinite comparison
--    let (inSet, outSet) = partition (b x) l
        ass' = ass ++ [inSet]
    in  pGroupPKernelByTerm b outSet ass'

pPartition :: (PKernel -> PKernel -> Bool) -> PKernel -> [PKernel] -> ([PKernel], [PKernel]) -> ([PKernel], [PKernel])
pPartition b x [] (inSet, outSet) = (inSet, outSet)
pPartition b x l@(a:as) (inSet, outSet) =
    if b x a
       then let inSet' = inSet ++ [a]
            in  pPartition b a as (inSet', outSet)
       else let outSet' = outSet ++ [a]
            in  pPartition b x as (inSet, outSet')

pGroupBy :: (a -> a -> Bool) -> [a] -> [[a]]
pGroupBy b l = pGroupByTerm b l []

pGroupByTerm :: (a -> a -> Bool) -> [a] -> [[a]] -> [[a]]
pGroupByTerm b [] ass = ass
pGroupByTerm b l@(x:xs) ass =
    let (as1, as2) = partition (b x) l
        ass' = if null as1 then ass else ass ++ [as1]
    in  pGroupByTerm b as2 ass'

pPackKernelFuncsIntoMTiles :: [PKernelFunc] -> [PMTile]
pPackKernelFuncsIntoMTiles l@(x:xs) =
    let l_mtile = mkMTile x
    in  pPackKernelFuncsIntoMTilesTerm xs [l_mtile]
        where pPackKernelFuncsIntoMTilesTerm [] ts = ts
              pPackKernelFuncsIntoMTilesTerm l@(x:xs) ts = 
                  case kfTileOp x of
                      PNULL -> pPackKernelFuncsIntoMTilesTerm xs ts
                      otherwise -> pPackKernelFuncsIntoMTilesTerm xs $ pushBackMTile x ts

mkMTile :: PKernelFunc -> PMTile
mkMTile l_kernel_func =
    let l_indices = kfTileIndex l_kernel_func
        l_sizes = kfTileSizes l_kernel_func
        l_term = mkMTileTerm l_indices l_sizes (ST [l_kernel_func])
    in  PMTile { mtSizes = l_sizes, mtKernelFuncs = [l_kernel_func], mtTerms = [l_term] }

mkMTileTerm :: [Int] -> [Int] -> PMTileItem -> PMTileTerm
mkMTileTerm l_indices l_sizes l_item =
    let l_tileTerm = PMTileTerm { mttIndex = l_indices,
                                  mttSizes = l_sizes,
                                  mttItem = l_item }
    in  l_tileTerm

pushBackMTile :: PKernelFunc -> [PMTile] -> [PMTile]
pushBackMTile l_kernel_func [] =
    let l_mtile = mkMTile l_kernel_func
    in  [l_mtile]
pushBackMTile l_kernel_func tL@(t:ts) =
    let l_indices = kfTileIndex l_kernel_func
        l_sizes = kfTileSizes l_kernel_func
        l_t_sizes = mtSizes t
        l_terms = mtTerms t
        -- lcp : longest common prefix
        l_lcp = pLongestCommonPrefix l_sizes l_t_sizes
        l_len_lcp = length l_lcp
        l_len_sizes = length l_sizes
        l_len_t_sizes = length l_t_sizes
    in  if l_len_lcp < l_len_sizes && l_len_lcp < l_len_t_sizes
           then t:(pushBackMTile l_kernel_func ts)
           else let t' = t { mtSizes = if length l_sizes > length l_t_sizes 
                                          then l_sizes 
                                          else l_t_sizes
                             ,
                             mtKernelFuncs = mtKernelFuncs t ++ [l_kernel_func],
                             mtTerms = pushBackMTileTerm l_indices l_sizes l_kernel_func l_terms }
                in  t':ts

pushBackMTileTerm :: [Int] -> [Int] -> PKernelFunc -> [PMTileTerm] -> [PMTileTerm]
pushBackMTileTerm l_indices l_sizes l_kernel_func [] = []
pushBackMTileTerm l_indices l_sizes l_kernel_func l_terms@(t:ts) =
    let l_t_indices = mttIndex t
        l_lcp = pLongestCommonPrefix l_indices l_t_indices
        l_len_lcp = length l_lcp
        l_len_t = length l_t_indices
        l_len_indices = length l_indices
    in  if l_len_lcp == l_len_t
           -- 'then' implies that length l_indices >= length l_t_indices
           then let l_mitem = pushBackMTileItem (drop l_len_lcp l_indices) (drop l_len_lcp l_sizes) l_kernel_func (mttItem t)
                    t' = t { mttItem = l_mitem }
                -- in  t':(pushBackMTileTerm l_indices l_sizes l_kernel_func ts)
                in  t':ts
           -- else: implies length l_indices < length l_t_indices 
           --       we have to split 
           else if l_len_lcp == l_len_indices || (l_len_lcp > 0 && null ts)
                   then let l_this_index = take l_len_lcp $ mttIndex t
                            l_this_sizes = take l_len_lcp $ mttSizes t
                            l_new_mterm_0 = PMTileTerm {
                                mttIndex = (drop l_len_lcp $ mttIndex t),
                                mttSizes = (drop l_len_lcp $ mttSizes t),
                                mttItem = mttItem t }
                            l_new_mterm_1 = PMTileTerm {
                                mttIndex = (drop l_len_lcp l_indices),
                                mttSizes = (drop l_len_lcp l_sizes),
                                mttItem = (ST [l_kernel_func]) }
                            t' = PMTileTerm {
                                mttIndex = l_this_index,
                                mttSizes = l_this_sizes,
                                mttItem = MT ([l_new_mterm_0] ++ [l_new_mterm_1]) }
                        in  t':(pushBackMTileTerm l_indices l_sizes l_kernel_func ts)
                   -- else : implies l_len_lcp < l_len_t && l_len_lcp < l_len_indices 
                   --        we have to split
                   else if (l_len_lcp == 0 && null ts)
                           then let t' = mkMTileTerm l_indices l_sizes (ST [l_kernel_func])
                                in  [t] ++ [t']
                           else t:(pushBackMTileTerm l_indices l_sizes l_kernel_func ts)

pushBackMTileItem :: [Int] -> [Int] -> PKernelFunc -> PMTileItem -> PMTileItem
pushBackMTileItem [] [] l_kernel_func (ST l_ks) = 
    ST (l_ks ++ [l_kernel_func])
pushBackMTileItem l_indices l_sizes l_kernel_func (ST l_ks) =
    let l_old_term = PMTileTerm {
            mttIndex = [],
            mttSizes = [],
            mttItem = ST l_ks }
        l_new_term = PMTileTerm {
            mttIndex = l_indices,
            mttSizes = l_sizes,
            mttItem = ST [l_kernel_func] }
    in  MT ([l_old_term] ++ [l_new_term])
pushBackMTileItem l_indices l_sizes l_kernel_func (MT l_ts) =
    let l_new_terms = pushBackMTileTerm l_indices l_sizes l_kernel_func l_ts
    in  MT l_new_terms

pSetTileOp :: TileOp -> PTile -> PTile
pSetTileOp l_op pTile = pTile { tOp = l_op }

pStripPrefixUnderScore :: String -> String
pStripPrefixUnderScore l_str =
    if isPrefixOf "__" l_str
        then drop 2 l_str
        else l_str

pStripSuffixUnderScore :: String -> String
pStripSuffixUnderScore l_str =
    if isSuffixOf "__" l_str
        then take (length l_str - 2) l_str
        else l_str

pFillGuardOrder :: Int -> PGuard -> PGuard
pFillGuardOrder n l_guard = l_guard { gOrder = n }

pFillTileOrder :: Int -> PTile -> PTile
pFillTileOrder n l_tile = l_tile { tOrder = n }

pFillPShapeName :: PName -> Int -> PShape -> PShape
pFillPShapeName l_name n l_pShape = l_pShape { shapeName = "__POCHOIR_Shape_" ++ l_name ++ "_" ++ show n ++ "__" }

pGetOverlapGuardFuncName :: PGuard -> String
pGetOverlapGuardFuncName l_pGuard = 
    let l_gComment = gComment l_pGuard
        l_gOrder = gOrder l_pGuard
        l_gfName = filter (not . isPrefixOf "!") l_gComment
        l_gfName' = "__POCHOIR_Guard_" ++ show l_gOrder ++ "__"
    in  l_gfName'

pGetOverlapGuardName :: PGuard -> String
pGetOverlapGuardName l_pGuard =
    let l_gComment = gComment l_pGuard
        l_gOrder = gOrder l_pGuard
        l_gfName = filter (not . isPrefixOf "!") l_gComment
        l_gfName' = "__POCHOIR_Guard_" ++ show l_gOrder ++ "__"
        l_gName = pStripSuffixUnderScore $ pStripPrefixUnderScore l_gfName' 
    in  l_gName

pGetAllIGuardTiles :: Int -> [String] -> [(PGuard, PTile)] -> [PTile] -> [(PGuard, [PTile])]
pGetAllIGuardTiles l_rank l_condStr [] l_tiles =
    let l_pGuard = PGuard { gName = "__" ++ (intercalate "_" l_condStr) ++ "__", gRank = l_rank, gFunc = emptyGuardFunc, gComment = l_condStr, gOrder = 0, gColor = emptyColor }
        l_tiles' = map (pSetTileOp PSERIAL) l_tiles
    in  [(l_pGuard, l_tiles')]
pGetAllIGuardTiles l_rank l_condStr l_iGTs@(i:is) l_tiles =
    let l_condStr' = l_condStr ++ [gName $ fst i]
        l_tiles' = l_tiles ++ [snd i]
        l_condStr'' = l_condStr ++ ["!" ++ (gName $ fst i)]
        l_tiles'' = l_tiles
    in  pGetAllIGuardTiles l_rank l_condStr' is l_tiles' ++ pGetAllIGuardTiles l_rank l_condStr'' is l_tiles''

pInvertIdxN :: Int -> Int -> [String ]
pInvertIdxN l_xIdx l_xLen =
    let l_condStr = if l_xIdx < l_xLen 
                       then replicate l_xIdx "!" ++ [""] ++ replicate (l_xLen - l_xIdx - 1) "!"
                       else []
    in  l_condStr

pGetExclusiveGuardTiles :: PMode -> Int -> [(PGuard, PTile)] -> [(PGuard, PTile)] -> [(PGuard, PTile)] -> [(PGuard, [PTile])]
pGetExclusiveGuardTiles l_mode l_rank [] l_iGTs l_tiGTs = 
    let l_condStr = if null l_tiGTs then [] else map ((++) "!" . gName . fst) l_tiGTs 
        l_pGuardTiles = pGetAllIGuardTiles l_rank l_condStr l_iGTs []
    in  l_pGuardTiles 
pGetExclusiveGuardTiles l_mode l_rank l_xGTs l_iGTs l_tiGTs =
    let l_condStr = if null l_tiGTs then map ((++) "!" . gName . fst) l_xGTs else map ((++) "!" . gName . fst) l_xGTs ++ map ((++) "!" . gName . fst) l_tiGTs 
        l_pGuardTiles = pGetAllIGuardTiles l_rank l_condStr l_iGTs []
        l_xLen = length l_xGTs
    in  l_pGuardTiles ++ pGetExclusiveGuardTilesTerm l_mode l_rank 0 (l_xLen-1) l_xGTs l_iGTs l_tiGTs

pGetExclusiveGuardTilesTerm :: PMode -> Int -> Int -> Int -> [(PGuard, PTile)] -> [(PGuard, PTile)] -> [(PGuard, PTile)] -> [(PGuard, [PTile])]
pGetExclusiveGuardTilesTerm l_mode l_rank l_xIdx l_xLen [] l_iGTs l_tiGTs = []
pGetExclusiveGuardTilesTerm l_mode l_rank l_xIdx l_xLen l_xGTs l_iGTs l_tiGTs =
    let l_tiCondStr = if null l_tiGTs then [] else map ((++) "!" . gName . fst) l_tiGTs
        l_xCondStr = zipWith (++) (pInvertIdxN l_xIdx l_xLen) (map (gName . fst) l_xGTs)
        l_x = l_xGTs !! l_xIdx
        l_condStr = l_xCondStr ++ l_tiCondStr
        l_pGuardTiles = if l_xIdx < l_xLen 
                           then pGetAllIGuardTiles l_rank l_condStr l_iGTs [snd l_x]
                           else []
        l_rest = if l_xIdx < l_xLen
                    then pGetExclusiveGuardTilesTerm l_mode l_rank (l_xIdx + 1) l_xLen l_xGTs l_iGTs l_tiGTs
                    else []
    in  l_pGuardTiles ++ l_rest

pGetInclusiveGuardTiles :: PMode -> [(PGuard, PTile)] -> [(PGuard, PTile)] -> [(PGuard, PTile)] -> [(PGuard, [PTile])]
pGetInclusiveGuardTiles _ _ _ [] = [] 
pGetInclusiveGuardTiles l_mode [] [] l_tiGTs =
    let l_tiGs = map fst l_tiGTs
        l_iTs = map (pSetTileOp PINCLUSIVE . snd) l_tiGTs
        l_tigStr = map gName l_tiGs
        l_rank = if null l_tiGs then 0 else (gRank . head) l_tiGs
        l_pIG = emptyGuard { gName = "__" ++ (intercalate "_" l_tigStr) ++ "__", gRank = l_rank, gComment = l_tigStr }
    in  [(l_pIG, l_iTs)]
pGetInclusiveGuardTiles l_mode [] l_iGTs l_tiGTs =
    let l_tiGs = map fst l_tiGTs
        l_iGs = map fst l_iGTs
        l_iTs = map (pSetTileOp PINCLUSIVE . snd) (l_iGTs ++ l_tiGTs)
        l_tigStr = map gName l_tiGs
        l_rank = if null l_tiGs 
                    then if null l_iGTs 
                            then 0 
                            else (gRank . head) l_iGs 
                    else (gRank . head) l_tiGs
        l_pIG = emptyGuard { gName = "__" ++ (intercalate "_" l_tigStr) ++ "__", gRank = l_rank, gComment = l_tigStr }
    in  [(l_pIG, l_iTs)]
pGetInclusiveGuardTiles l_mode l_xGTs l_iGTs l_tiGTs =
    let l_tiGs = map fst l_tiGTs
        l_iGs = map fst l_iGTs
        l_xTs = map (pSetTileOp PEXCLUSIVE . snd) l_xGTs
        l_iTs = map (pSetTileOp PINCLUSIVE . snd) (l_iGTs ++ l_tiGTs)
        l_tigStr = map gName l_tiGs
        l_rank = if null l_tiGs 
                    then if null l_iGTs 
                            then 0 
                            else (gRank . head) l_iGs 
                    else (gRank . head) l_tiGs
        l_pIG = emptyGuard { gName = "__" ++ (intercalate "_" l_tigStr) ++ "__", gRank = l_rank, gComment = l_tigStr }
    in  [(l_pIG, l_xTs ++ l_iTs)]

pAddUnderScore :: String -> String
pAddUnderScore str =
    if isPrefixOf "!" str
        then "!__" ++ drop 1 str ++ "__"
        else "__" ++ str ++ "__"

pSubstitute :: Eq a => [a] -> [a] -> [a] -> [a]
pSubstitute _ _ [] = []
pSubstitute from to xs@(a:as) =
    if isPrefixOf from xs
        then to ++ (pSubstitute from to $ drop (length from) xs)
        else a:pSubstitute from to as

getPShape :: ParserState -> String -> PShape
getPShape l_state l_shape =
    case Map.lookup l_shape $ pShape l_state of
        Nothing -> emptyShape
        Just l_pShape -> l_pShape

getPStencil :: String -> ParserState -> PStencil -> PStencil
getPStencil l_id l_state l_oldStencil =
    case Map.lookup l_id $ pStencil l_state of
        Nothing -> l_oldStencil
        Just l_stencil -> l_stencil

-- generate new kernel parameters (t, i0, i1, ...) according to the number
-- of old kernel params
getKernelParams :: Int -> [String]
getKernelParams 0 = []
getKernelParams n = ["t"] ++ getSpaceKernelParams (n-1)
    where getSpaceKernelParams 0 = []
          getSpaceKernelParams k = let i' = "i" ++ show (k-1)
                                   in  [i'] ++ getSpaceKernelParams (k-1)

getIter :: PArray -> Expr -> PRWMode -> [Iter]
getIter arrayInUse (PVAR q v dL) l_rw =
    let iterName = "iter"
    in  [(iterName, arrayInUse, dL, l_rw)]
getIter arrayInUse _ _ = []

getBaseIter :: [PName] -> [Iter] -> [Iter]
getBaseIter _ []  = []
getBaseIter l_kernelParams iL@(i:is) = union (getBaseIterItem i) (getBaseIter l_kernelParams is)
    where getBaseIterItem (name, array, dims, rw) = 
            let dims' = getBaseDimExpr (tail l_kernelParams) dims
            in  [("baseIter_", array, dims', rw)]

getPointer :: [PName] -> PArray -> Expr -> PRWMode -> [Iter]
getPointer kernelParams arrayInUse (PVAR q v dL) l_rw =
    let iterName = "pt_" ++ aName arrayInUse ++ "_"
        dL' = transDimExpr kernelParams dL
    in  [(iterName, arrayInUse, dL', l_rw)]
getPointer _ arrayInUse _ _ = []

transIterN :: Int -> [Iter] -> [Iter]
transIterN _ [] = []
transIterN n ((name, array, dim, rw):is) = (name ++ show n, array, dim, rw) : (transIterN (n+1) is)

pFillIters :: [Iter] -> PKernelFunc -> PKernelFunc
pFillIters l_iters l_kernel_func = l_kernel_func { kfIter = l_iters}

pTileLength :: PTile -> Int
pTileLength l_tile = (head . tSizes) l_tile
{-
    let l_tile_kernel = tKernel l_tile
    in  case l_tile_kernel of
            SK _ -> 1
            LK l_tKs -> length l_tKs
-}

transArrayMap :: [PArray] -> Map.Map PName PArray
transArrayMap aL = Map.fromList $ transAssocList aL
    where transAssocList [] = []
          transAssocList aL@(a:as) = [(aName a, a)] ++ (transAssocList as)

getBaseDimExpr :: [PName] -> [DimExpr] -> [DimExpr]
getBaseDimExpr kL [] = []
getBaseDimExpr kL dL@(d:ds) = [d] ++ map transKernelParamsToDimExprItem kL
    where transKernelParamsToDimExprItem k = DimVAR k

transDimExpr :: [PName] -> [DimExpr] -> [DimExpr]
transDimExpr kL [] = []
transDimExpr kL dL@(d:ds) = [d] ++ (transSpaceDimExpr ds)
    where transSpaceDimExpr [] = []
          transSpaceDimExpr (d:ds) = (transSpaceDimExprItem d) ++ (transSpaceDimExpr ds)
          transSpaceDimExprItem (DimParen e) =
                        transSpaceDimExprItem e
          transSpaceDimExprItem (DimVAR v) =
                        if elem v kL then [(DimVAR v)]
                                     else []
          transSpaceDimExprItem (DimDuo bop e1 e2) =
                        transSpaceDimExprItem e1 ++ transSpaceDimExprItem e2
          transSpaceDimExprItem (DimINT n) = []

pFirst (a, _, _) = a
pSecond (_, b, _) = b
pThird (_, _, c) = c

pIterName (a, _, _, _) = a
pIterArray (_, a, _, _) = a
pIterDims (_, _, a, _) = a
pIterRWMode (_, _, _, a) = a

transInterior :: [PName] -> Expr -> Expr
transInterior l_arrayInUse (PVAR q v dL) =
    if elem v l_arrayInUse == True then PVAR q (v ++ ".interior") dL
                                   else PVAR q v dL
transInterior l_arrayInUse e = e

getArrayName :: [PArray] -> [PName]
getArrayName [] = []
getArrayName (a:as) = (aName a) : (getArrayName as)

pGetReadIters :: [Iter] -> [Iter]
pGetReadIters [] = []
pGetReadIters iL@(i:is) = 
    if isReadIter i == True then union [i] (pGetReadIters is) else pGetReadIters is
       where isReadIter (iName, iArray, iDims, iRW) =
               if iRW == PRead then True else False

pGetWriteIters :: [Iter] -> [Iter]
pGetWriteIters [] = []
pGetWriteIters iL@(i:is) = 
    if isWriteIter i == True then union [i] (pGetWriteIters is) else pGetWriteIters is
       where isWriteIter (iName, iArray, iDims, iRW) =
               if iRW == PWrite then True else False

pGetMinIters :: [Iter] -> [Iter]
pGetMinIters iL@(i:is) = map (minimumBy cmpDim) $ groupIterByArray $ sortIterByArray iL
    where cmpDim (iName, iArray, iDims, iRW) (jName, jArray, jDims, jRW) = compare iDims jDims 

groupIterByArray :: [Iter] -> [[Iter]]
groupIterByArray iL@(i:is) = pGroupBy eqArray iL
    where eqArray (iName, iArray, iDims, iRW) (jName, jArray, jDims, jRW) = aName iArray == aName jArray 

sortIterByArray :: [Iter] -> [Iter]
sortIterByArray iL@(i:is) = sortBy cmpArray iL
    where cmpArray (iName, iArray, iDims, iRW) (jName, jArray, jDims, jRW) = compare (aName iArray) (aName jArray)

mergePShapes :: PShape -> PShape -> PShape
mergePShapes a b = PShape { shapeName = shapeName a ++ shapeName b,
                            shapeRank = max (shapeRank a) (shapeRank b),
                            shapeLen = length $ union (shape a) (shape b),
                            shapeToggle = max (shapeToggle a) (shapeToggle b),
                            shapeSlopes = zipWith max (shapeSlopes a) (shapeSlopes b),
                            shapeTimeShift = min (shapeTimeShift a) (shapeTimeShift b),
                            shape = union (shape a) (shape b),
                            shapeComment = "" }

pSysShape :: PShape -> PShape
pSysShape a = a { shapeName = "__" ++ shapeName a ++ "__" }

pSys :: String -> String
pSys a = "__" ++ a ++ "__"

checkValidPArray :: String -> ParserState -> (Bool, PArray)
checkValidPArray l_array l_state =
    case Map.lookup l_array $ pArray l_state of
         Nothing -> (False, emptyPArray)
         Just l_pArray -> (True, l_pArray)

fillToggleInPArray :: PArray -> Int -> PArray
fillToggleInPArray l_pArray l_toggle = l_pArray { aToggle = l_toggle }

-- convert a binary number (in String) into a decimal number (in Int)
pBinToDec :: String -> Int
pBinToDec bValue = pBinToDecTerm 0 (length bValue) bValue
    where pBinToDecTerm idx n bValue =
            if (idx >= n) 
               then 0
               else if bValue !! idx == '1'
                       then (shift 1 $ n - 1 - idx) + pBinToDecTerm (idx+1) n bValue
                       else pBinToDecTerm (idx+1) n bValue

rename :: String -> String -> String
rename pSuffix fname = name ++ pSuffix ++ ".cpp"
    where (name, suffix) = break ('.' ==) fname

pInsMod :: PName -> Int -> Int -> String
pInsMod l_param l_index l_size =
    (mkParen $ l_param ++ " % " ++ show l_size) ++ " == " ++ show l_index

mkInput :: String -> String
mkInput a = "_" ++ a

mkOutput :: String -> String
mkOutput a = a ++ "_"

mkLocal :: String -> String
mkLocal a = "l_" ++ a

mkParen :: String -> String
mkParen a = "(" ++ a ++ ")"

mkComment :: String -> String
mkComment a = "/* " ++ a ++ " */ "

mkPointer :: String -> String
mkPointer a = a ++ " * "

derefPointer :: String -> String
derefPointer a = "*" ++ a

mkStatic :: String -> String
mkStatic a = "static " ++ a
-- mkStatic a = a

pDeclVar :: String -> String -> String -> String
pDeclVar l_type l_name l_init =
    if null l_init 
       then l_type ++ l_name ++ ";"
       else l_type ++ l_name ++ " = " ++ l_init ++ ";"

pAssign :: String -> String -> String
pAssign l_name l_value =
    l_name ++ " = " ++ l_value ++ ";"

zipInsert :: String -> String -> String -> String -> String
zipInsert op delim a b = a ++ op ++ b ++ delim

getTagFromMode :: PMode -> String
getTagFromMode l_mode = 
    case l_mode of
        PAllCondTileMacroOverlap -> "Macro"
        PAllCondTileCPointerOverlap -> "CPointer"
        PAllCondTilePointerOverlap -> "Pointer"
        PAllCondTileOptPointerOverlap -> "Opt_Pointer"
        PUnrollTimeTileMacroOverlap -> "Macro"
        PUnrollTimeTileCPointerOverlap -> "CPointer"
        PUnrollTimeTilePointerOverlap -> "Pointer"
        PUnrollTimeTileOptPointerOverlap -> "Opt_Pointer"
        otherwise -> "Unknown"

foldTileOp :: TileOp -> TileOp -> TileOp
foldTileOp PNULL PNULL = PNULL
foldTileOp _ _ = PSERIAL

pIsTileOpNull :: TileOp -> Bool
pIsTileOpNull PNULL = True
pIsTileOpNull _ = False
