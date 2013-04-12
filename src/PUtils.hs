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

module PUtils where

import Text.ParserCombinators.Parsec
import Control.Monad
import Data.List
import PData
import qualified Data.Map as Map

updatePTile :: PTile -> ParserState -> ParserState
updatePTile l_tile parserState =
    parserState { pTile = Map.insert (tName l_tile) l_tile $ pTile parserState }

updatePTileOrder :: Int -> ParserState -> ParserState
updatePTileOrder l_tile_order parserState =
    parserState { pTileOrder = l_tile_order }

updatePColorNum :: Int -> ParserState -> ParserState
updatePColorNum l_color_num parserState =
    parserState { pColorNum = l_color_num }

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

updatePMacro :: (PName, PValue) -> ParserState -> ParserState
updatePMacro (l_name, l_value) parserState =
    parserState { pMacro = Map.insert l_name l_value (pMacro parserState) }

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

updateStencilRegStaggerKernel :: String -> [(PGuard, [PKernel])] -> ParserState -> ParserState
updateStencilRegStaggerKernel l_id l_regStaggerKernel parserState =
    let f k x =
            if sName x == l_id
                then Just $ x { sRegStaggerKernel = l_regStaggerKernel }
                else Nothing
    in  parserState { pStencil = Map.updateWithKey f l_id $ pStencil parserState }

updateStencilRegTileKernel :: String -> [(PGuard, PTile)] -> ParserState -> ParserState
updateStencilRegTileKernel l_id l_regTileKernel parserState =
    let f k x =
            if sName x == l_id
                then Just $ x { sRegTileKernel = l_regTileKernel }
                else Nothing
    in  parserState { pStencil = Map.updateWithKey f l_id $ pStencil parserState }

updateStencilRegInclusiveTileKernel :: String -> [(PGuard, PTile)] -> ParserState -> ParserState
updateStencilRegInclusiveTileKernel l_id l_regInclusiveTileKernel parserState =
    let f k x =
            if sName x == l_id
                then Just $ x { sRegInclusiveTileKernel = l_regInclusiveTileKernel }
                else Nothing
    in  parserState { pStencil = Map.updateWithKey f l_id $ pStencil parserState }

updateStencilRegTinyInclusiveTileKernel :: String -> [(PGuard, PTile)] -> ParserState -> ParserState
updateStencilRegTinyInclusiveTileKernel l_id l_regTinyInclusiveTileKernel parserState =
    let f k x =
            if sName x == l_id
                then Just $ x { sRegTinyInclusiveTileKernel = l_regTinyInclusiveTileKernel }
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
    let l_t = map pHead l_shapes
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
    let l_t = map pHead l_shapes
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
                                    Just l_kernel -> l_kernel
                l_tileKernel = SK l_kernel
            in  PTile { tName = l_tileName, tRank = kRank l_kernel, tSize = [1], tKernel = l_tileKernel, tComment = "", tOp = PNOP, tOrigGuard = emptyGuard, tOrder = 0 }
        Just l_pTile -> l_pTile { tComment = cKnown "Tile" }

getTileKernels :: PTile -> [PKernel]
getTileKernels l_tile = getKernelsTile [] 0 l_tile $ tKernel l_tile

getKernelsTile :: [Int] -> Int -> PTile -> PTileKernel -> [PKernel]
getKernelsTile l_indices l_index l_tile (SK l_kernel) =
    let ll_indices = l_indices ++ [l_index]
        ll_kernel_func = kFunc l_kernel
        ll_indices' = if null ll_indices then [] else init ll_indices
        ll_indices'' = if null ll_indices' then [0] else ll_indices'
        l_gfunc = gFunc $ tOrigGuard l_tile
        l_tile_op = tOp l_tile
        l_tile_order = tOrder l_tile
        ll_kernel = l_kernel { kIndex = ll_indices'', kFunc = ll_kernel_func { kfGuardFunc = l_gfunc, kfTileOp = l_tile_op, kfTileOrder = l_tile_order } }
    in  [ll_kernel]
getKernelsTile l_indices l_index l_tile (LK l_tKs@(t:ts)) =
    let ll_indices = l_indices ++ [l_index]
    in  getKernelsTile ll_indices 0 l_tile t ++ 
        getKernelsTile l_indices (l_index + 1) l_tile (LK ts)
getKernelsTile l_indices l_index l_tile (LK []) = []

getTileSizes :: [[Int]] -> [Int]
getTileSizes l_tile_indices = map ((+) 1 . maximum) $ transpose l_tile_indices

eqTileOpPKernelFunc :: PKernelFunc -> PKernelFunc -> Bool
eqTileOpPKernelFunc a b = (kfTileOp a) == (kfTileOp b)

eqTPKernel :: PKernel -> PKernel -> Bool
eqTPKernel a b = (pHead $ kIndex a) == (pHead $ kIndex b)

eqTGroupPKernel :: [PKernel] -> [PKernel] -> Bool
eqTGroupPKernel [] [] = True
eqTGroupPKernel a b = (pHead . kIndex . pHead) a == (pHead . kIndex . pHead) b

eqTPTile :: [Int] -> [Int] -> Bool
eqTPTile [] [] = True
eqTPTile [] b = False
eqTPTile a [] = False
eqTPTile a b = pHead a == pHead b

eqTileIndex :: [Int] -> [Int] -> Bool
eqTileIndex [] [] = True
eqTileIndex _ [] = True
eqTileIndex [] _ = True
eqTileIndex aL@(a:as) bL@(b:bs) = 
    if a == b
       then eqTileIndex as bs
       else False

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
        l_tOrder = kfTileOrder $ kFunc k
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

pPartition :: (PKernel -> PKernel -> Bool) -> PKernel -> [PKernel] -> ([PKernel], [PKernel]) -> ([PKernel], [PKernel])
pPartition b x [] (inSet, outSet) = (inSet, outSet)
pPartition b x l@(a:as) (inSet, outSet) =
    if b x a
       then let inSet' = inSet ++ [a]
            in  pPartition b a as (inSet', outSet)
       else let outSet' = outSet ++ [a]
            in  pPartition b x as (inSet, outSet')

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
    if ((length . kIndex . last . last) x <= (length . kIndex . pHead . pHead) y 
        && length x <= length y)
        then let x' = pMergeForward b x y
             in  pGroupPKernelByMergeTerm b (x':zs)
        else if ((length . kIndex . last . last) x > (length . kIndex . pHead . pHead) y
               && length x > length y)
                then let x' = pMergeBackward b x y []
                     in  pGroupPKernelByMergeTerm b (x':zs)
                else x ++ pGroupPKernelByMergeTerm b (y:zs)

pMergeForward :: (PKernel -> PKernel -> Bool) -> [[PKernel]] -> [[PKernel]] -> [[PKernel]]
pMergeForward b [] y = y
pMergeForward b x@(a:as) y = 
    let could_merge = foldr (||) False $ 
                        map (b (last a)) (map pHead y)
    in  if could_merge
           then pMergeForward b as $ pMergeForwardTerm b a y
           else a:(pMergeForward b as y)

pMergeForwardTerm :: (PKernel -> PKernel -> Bool) -> [PKernel] -> [[PKernel]] -> [[PKernel]]
pMergeForwardTerm b x [] = []
pMergeForwardTerm b x yL@(y:ys) =
    if b (last x) (pHead y)
       then let x' = map (pFillKIndex $ (kIndex . pHead) y) x
            in  (x' ++ y):(pMergeForwardTerm b x ys) 
       else y:(pMergeForwardTerm b x ys)

pMergeBackward :: (PKernel -> PKernel -> Bool) -> [[PKernel]] -> [[PKernel]] -> [[PKernel]] -> [[PKernel]]
pMergeBackward b x [] l = x ++ l
pMergeBackward b x y@(a:as) l =
    let could_merge = foldr (||) False $
                        map (flip b (pHead a)) (map last x)
    in  if could_merge 
           then pMergeBackward b (pMergeBackwardTerm b x a) as l
           else pMergeBackward b x as $ l ++ [a]

pMergeBackwardTerm :: (PKernel -> PKernel -> Bool) -> [[PKernel]] -> [PKernel] -> [[PKernel]]
pMergeBackwardTerm b [] y = []
pMergeBackwardTerm b xL@(x:xs) y =
    if b (last x) (pHead y)
       then let y' = map (pFillKIndex $ (kIndex . last) x) y
            in  (x ++ y'):(pMergeBackwardTerm b xs y)
       else x:(pMergeBackwardTerm b xs y)

pGroupPKernelBy :: (PKernel -> PKernel -> Bool) -> [PKernel] -> [[PKernel]]
pGroupPKernelBy b [] = []
pGroupPKernelBy b l = 
    let n = (kfTileOrder . kFunc . pHead) l
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

pGroupBy :: (a -> a -> Bool) -> [a] -> [[a]]
pGroupBy b l = pGroupByTerm b l []

pGroupByTerm :: (a -> a -> Bool) -> [a] -> [[a]] -> [[a]]
pGroupByTerm b [] ass = ass
pGroupByTerm b l@(x:xs) ass =
    let (as1, as2) = partition (b x) l
        ass' = if null as1 then ass else ass ++ [as1]
    in  pGroupByTerm b as2 ass'

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

pFillGuardOrder :: Int -> (PGuard, [PTile]) -> (PGuard, [PTile])
pFillGuardOrder n (l_pGuard, l_pTiles) = (l_pGuard { gOrder = n }, l_pTiles)

pFillPShapeName :: PName -> Int -> PShape -> PShape
pFillPShapeName l_name n l_pShape = l_pShape { shapeName = "__POCHOIR_Shape_" ++ l_name ++ "_" ++ show n ++ "__" }

pGetGuardFuncName :: PGuard -> String
pGetGuardFuncName l_pGuard = 
    let l_gComment = gComment l_pGuard
        l_gOrder = gOrder l_pGuard
        l_gfName = filter (not . isPrefixOf "!") l_gComment
        l_gfName' = "__POCHOIR_Guard_" ++ show l_gOrder ++ "__"
    in  l_gfName'

pGetGuardName :: PGuard -> String
pGetGuardName l_pGuard =
    let l_gComment = gComment l_pGuard
        l_gOrder = gOrder l_pGuard
        l_gfName = filter (not . isPrefixOf "!") l_gComment
        l_gfName' = "__POCHOIR_Guard_" ++ show l_gOrder ++ "__"
        l_gName = pStripSuffixUnderScore $ pStripPrefixUnderScore l_gfName' 
    in  l_gName

pGetAllIGuardTiles :: Int -> [String] -> [(PGuard, PTile)] -> [PTile] -> [(PGuard, [PTile])]
pGetAllIGuardTiles l_rank l_condStr [] l_tiles =
    let l_pGuard = PGuard { gName = "__" ++ (intercalate "_" l_condStr) ++ "__", gRank = l_rank, gFunc = emptyGuardFunc, gComment = l_condStr, gOrder = 0, gTileOrder = 0 }
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
    in  l_pGuardTiles ++ pGetExclusiveGuardTilesTerm l_mode l_rank 0 l_xLen l_xGTs l_iGTs l_tiGTs

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
        l_rank = if null l_tiGs then 0 else (gRank . pHead) l_tiGs
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
                            else (gRank . pHead) l_iGs 
                    else (gRank . pHead) l_tiGs
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
                            else (gRank . pHead) l_iGs 
                    else (gRank . pHead) l_tiGs
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
pTileLength l_tile =
    let l_tile_kernel = tKernel l_tile
    in  case l_tile_kernel of
            SK _ -> 1
            LK l_tKs -> length l_tKs

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

pHead :: [a] -> a
pHead = head

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

rename :: String -> String -> String
rename pSuffix fname = name ++ pSuffix ++ ".cpp"
    where (name, suffix) = break ('.' ==) fname

mkLocalFile :: String -> String 
mkLocalFile a = "./" ++ a

mkQuote :: String -> String
mkQuote a = "\"" ++ a ++ "\""
