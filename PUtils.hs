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

updatePKernel :: PKernel -> ParserState -> ParserState
updatePKernel l_kernel parserState =
    parserState { pKernel = Map.insert (kName l_kernel) l_kernel $ pKernel parserState }

updatePGuard :: PGuard -> ParserState -> ParserState
updatePGuard l_guard parserState =
    parserState { pGuard = Map.insert (gName l_guard) l_guard $ pGuard parserState }

updateObase :: PMode -> ParserState -> ParserState
updateObase mode parserState = parserState { pMode = mode }

updatePState :: PState -> ParserState -> ParserState
updatePState pstate parserState
    | (pstate == PochoirBegin && pState parserState == Unrelated) = parserState { pState = PochoirBegin }
    | (pstate == PochoirEnd && pState parserState == PochoirBegin) = parserState { pState = Unrelated }
    | ((pstate == PochoirMacro || pstate == PochoirDeclArray || pstate == PochoirDeclRange)
            && pState parserState == Unrelated) = parserState { pState = pstate }
    | pstate == Unrelated = parserState { pState = Unrelated }
    | otherwise = parserState { pState = PochoirError }

updatePMacro :: (PName, PValue) -> ParserState -> ParserState
updatePMacro (l_name, l_value) parserState =
    parserState { pMacro = Map.insert l_name l_value (pMacro parserState) }

updatePShape :: (PName, Int, PValue, Int, [Int], Int, [[Int]]) -> ParserState -> ParserState
updatePShape (l_name, l_rank, l_len, l_toggle, l_slopes, l_timeShift, l_shape) parserState =
    let l_pShape = PShape {shapeName = l_name, shapeRank = l_rank, shapeLen = l_len, shapeToggle = l_toggle, shapeSlopes = l_slopes, shapeTimeShift = l_timeShift, shape = l_shape}
    in parserState { pShape = Map.insert l_name l_pShape (pShape parserState) }

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

updateStencilArray :: String -> PArray -> ParserState -> ParserState
updateStencilArray l_id l_pArray parserState =
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

updateStencilUnroll :: String -> Int -> ParserState -> ParserState
updateStencilUnroll l_id l_unroll parserState =
    let f k x =
            if sName x == l_id
                then Just $ x { sUnroll = l_unroll }
                else Nothing
    in  parserState { pStencil = Map.updateWithKey f l_id $ pStencil parserState }

updateStencilToggle :: String -> Int -> ParserState -> ParserState
updateStencilToggle l_id l_toggle parserState =
    let f k x =
            if sName x == l_id
                then Just $ x { sToggle = l_toggle }
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
    
getValidKernel :: ParserState -> String -> (Bool, PKernel)
getValidKernel l_state l_kernel =
    case Map.lookup l_kernel $ pKernel l_state of
        Nothing -> (False, emptyKernel) 
        Just l_pKernel -> (True, l_pKernel)

getValidGuard :: ParserState -> String -> (Bool, PGuard)
getValidGuard l_state l_guard =
    case Map.lookup l_guard $ pGuard l_state of
        Nothing -> (False, emptyGuard)
        Just l_pGuard -> (True, l_pGuard)

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
groupIterByArray iL@(i:is) = groupBy eqArray iL
    where eqArray (iName, iArray, iDims, iRW) (jName, jArray, jDims, jRW) = aName iArray == aName jArray 

sortIterByArray :: [Iter] -> [Iter]
sortIterByArray iL@(i:is) = sortBy cmpArray iL
    where cmpArray (iName, iArray, iDims, iRW) (jName, jArray, jDims, jRW) = compare (aName iArray) (aName jArray)

