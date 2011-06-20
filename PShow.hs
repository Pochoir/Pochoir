{-
 ----------------------------------------------------------------------------------
 -  Copyright (C) 2010-2011  Massachusetts Institute of Technology
 -  Copyright (C) 2010-2011  Yuan Tang <yuantang@csail.mit.edu>
 - 		                     Charles E. Leiserson <cel@mit.edu>
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
module PShow where

import Text.ParserCombinators.Parsec
import Control.Monad
import qualified Data.Map as Map
import Data.List

import PData
import PUtils

simplifyDimExpr :: DimExpr -> DimExpr
simplifyDimExpr de = 
    let newDimExpr = simplifyDimExprItem de
    in  if newDimExpr == de then newDimExpr
                            else simplifyDimExpr newDimExpr

simplifyDimExprItem :: DimExpr -> DimExpr
simplifyDimExprItem (DimVAR v) = DimVAR v
simplifyDimExprItem (DimINT n) = DimINT n
simplifyDimExprItem (DimDuo "-" (DimINT m) (DimINT n)) = DimINT (m-n)
simplifyDimExprItem (DimDuo "+" (DimINT m) (DimINT n)) = DimINT (m+n)
simplifyDimExprItem (DimDuo "*" (DimINT m) (DimINT n)) = DimINT (m*n)
simplifyDimExprItem (DimDuo bop e1 e2)
    | bop == "+" && e1 == DimINT 0 = e2
    | bop == "+" && e2 == DimINT 0 = e1
    | bop == "-" && e2 == DimINT 0 = e1
    | bop == "*" && e1 == DimINT 0 = DimINT 0
    | bop == "*" && e2 == DimINT 0 = DimINT 0
    | bop == "*" && e1 == DimINT 1 = e2
    | bop == "*" && e2 == DimINT 1 = e1
    | otherwise = DimDuo bop (simplifyDimExprItem e1) (simplifyDimExprItem e2)
simplifyDimExprItem (DimParen e) = DimParen (simplifyDimExprItem e)

getFromStmts :: (PArray -> Expr -> PRWMode -> [Iter]) -> PRWMode -> Map.Map PName PArray -> [Stmt] -> [Iter]
getFromStmts l_action _ _ [] = []
getFromStmts l_action l_rw l_arrayMap l_stmts@(a:as) = 
    let i1 = getFromStmt a 
        i2 = getFromStmts l_action l_rw l_arrayMap as 
    in  union i1 i2
    where getFromStmt (BRACES stmts) = getFromStmts l_action l_rw l_arrayMap stmts 
          getFromStmt (EXPR e) = getFromExpr l_rw e
          getFromStmt (DEXPR qs t es) = concatMap (getFromExpr l_rw) es
          getFromStmt (IF e s1 s2) = 
              let iter1 = getFromExpr l_rw e 
                  iter2 = getFromStmt s1 
                  iter3 = getFromStmt s2 
              in  union iter1 (union iter2 iter3)
          getFromStmt (SWITCH e stmts) = 
              let iter1 = getFromExpr l_rw e 
                  iter2 = getFromStmts l_action l_rw l_arrayMap stmts
              in  union iter1 iter2
          getFromStmt (CASE v stmts) = getFromStmts l_action l_rw l_arrayMap stmts
          getFromStmt (DEFAULT stmts) = getFromStmts l_action l_rw l_arrayMap stmts
          getFromStmt (NOP) = []
          getFromStmt (BREAK) = []
          getFromStmt (DO e stmts) = 
              let iter1 = getFromExpr l_rw e 
                  iter2 = getFromStmts l_action l_rw l_arrayMap stmts
              in  union iter1 iter2
          getFromStmt (WHILE e stmts) = 
              let iter1 = getFromExpr l_rw e 
                  iter2 = getFromStmts l_action l_rw l_arrayMap stmts
              in  union iter1 iter2
          getFromStmt (FOR sL s) = 
              let iter1 = getFromStmt s 
                  iter2 = concat $ map (getFromStmts l_action l_rw l_arrayMap) sL
              in  union iter1 iter2
          getFromStmt (CONT) = []
          getFromStmt (RET e) = getFromExpr l_rw e
          getFromStmt (RETURN) = []
          getFromStmt (UNKNOWN s) = []
          getFromExpr _ (VAR q v) = []
          getFromExpr l_rw (PVAR q v dL) = 
              case Map.lookup v l_arrayMap of
                   Nothing -> []
                   Just arrayInUse -> l_action arrayInUse (PVAR q v dL) l_rw
          getFromExpr _ (BVAR v dim) = []
          getFromExpr l_rw (BExprVAR v e) = getFromExpr l_rw e
          getFromExpr l_rw (SVAR t e c f) = getFromExpr l_rw e
          getFromExpr l_rw (PSVAR t e c f) = getFromExpr l_rw e
          getFromExpr l_rw (Uno uop e) = getFromExpr l_rw e
          getFromExpr l_rw (PostUno uop e) = getFromExpr l_rw e
          getFromExpr l_rw (Duo bop e1 e2) = 
            if bop == "=" 
                then let iter1 = getFromExpr PWrite e1
                         iter2 = getFromExpr PRead e2
                     in  (union iter1 iter2)
                else let iter1 = getFromExpr l_rw e1
                         iter2 = getFromExpr l_rw e2
                     in  (union iter1 iter2)
          getFromExpr l_rw (PARENS e) = getFromExpr l_rw e
          getFromExpr _ _ = []

transStmts :: [Stmt] -> (Expr -> Expr) -> [Stmt]
transStmts [] _ = []
transStmts l_stmts@(a:as) l_action = transStmt a : transStmts as l_action
    where transStmt (BRACES stmts) = BRACES $ transStmts stmts l_action
          transStmt (EXPR e) = EXPR $ transExpr e
          transStmt (DEXPR qs t es) = DEXPR qs t $ map transExpr es 
          transStmt (IF e s1 s2) = IF (transExpr e) 
                                              (transStmt s1) 
                                              (transStmt s2) 
          transStmt (SWITCH e stmts) = SWITCH (transExpr e) 
                                                  (transStmts stmts l_action)
          transStmt (CASE v stmts) = CASE v $ transStmts stmts l_action
          transStmt (DEFAULT stmts) = DEFAULT $ transStmts stmts l_action
          transStmt NOP =  NOP
          transStmt BREAK =  BREAK
          transStmt (DO e stmts) = DO (transExpr e) (transStmts stmts l_action)
          transStmt (WHILE e stmts) = WHILE (transExpr e)
                                                (transStmts stmts l_action)
          transStmt (FOR sL s) = FOR (map (flip transStmts l_action) sL)
                                         (transStmt s) 
          transStmt (CONT) = CONT
          transStmt (RETURN) = RETURN
          transStmt (RET e) = RET (transExpr e)
          transStmt (UNKNOWN s) = UNKNOWN s 
          transExpr (VAR q v) =  VAR q v
          -- if it's in the form of BVAR, then the user must have already done some
          -- manual transformation on its source, we just leave them untouched!
          transExpr (PVAR q v dL) = l_action (PVAR q v dL)
          transExpr (BVAR v dim) = BVAR v dim
          transExpr (BExprVAR v e) = BExprVAR v $ transExpr e
          transExpr (SVAR t e c f) = SVAR t (transExpr e) c f
          transExpr (PSVAR t e c f) = PSVAR t (transExpr e) c f
          transExpr (Uno uop e) = Uno uop $ transExpr e
          transExpr (PostUno uop e) = PostUno uop $ transExpr e
          transExpr (Duo bop e1 e2) = Duo bop (transExpr e1) (transExpr e2)
          transExpr (PARENS e) = PARENS $ transExpr e
          transExpr (INT n) = (INT n)
          transExpr (FLOAT f) = (FLOAT f)
          transExpr (BOOL b) = (BOOL b)

pIterLookup :: (String, [DimExpr]) -> [Iter] -> Maybe String
pIterLookup (v, dL) [] = Nothing
pIterLookup (v, dL) ((iterName, arrayInUse, dim, rw):is) 
    | v == aName arrayInUse && dL == dim = Just iterName
    | otherwise = pIterLookup (v, dL) is

pShowShadowArrayInUse :: [PArray] -> String
pShowShadowArrayInUse [] = ""
pShowShadowArrayInUse aL@(a:as) =
    pShowShadowArrayItem a ++ pShowShadowArrayInUse as
    where pShowShadowArrayItem a = 
            let l_type = aType a 
                l_rank = aRank a
                l_toggle = aToggle a
                l_name = aName a
                pShowShadowHeader (l_type, l_rank, l_toggle) = "interior_shadow<" ++
                    show l_type ++ ", " ++ show l_rank ++ ", " ++ show l_toggle ++ "> "
            in breakline ++ pShowShadowHeader (l_type, l_rank, l_toggle) ++ 
                l_name ++ "_shadow(" ++ l_name ++ ");" ++
                breakline ++ pShowShadowHeader (l_type, l_rank, l_toggle) ++
                l_name ++ "(" ++ l_name ++ "_shadow);" ++ breakline

pDefMacroArrayInUse :: PName -> [PArray] -> [PName] -> String
pDefMacroArrayInUse _ [] _ = ""
pDefMacroArrayInUse l_macro (a:as) pL = pDefMacroShadowItem l_macro a pL ++ pDefMacroArrayInUse l_macro as pL
    where pDefMacroShadowItem l_macro a pL = 
            let l_arrayName = aName a
                l_arrayMacroName = l_arrayName ++ "." ++ l_macro
            in  "#define " ++ pShowArrayTerm l_arrayName pL ++ " " ++
                pShowArrayTerm l_arrayMacroName pL ++ breakline

pShowArrayTerm :: PName -> [PName] -> String
pShowArrayTerm a pL = a ++ "(" ++ pShowListIdentifiers pL ++ ")"

pUndefMacroArrayInUse :: [PArray] -> [PName] -> String
pUndefMacroArrayInUse [] _ = ""
pUndefMacroArrayInUse (a:as) pL = pUndefMacroShadowItem a pL ++ pUndefMacroArrayInUse as pL
    where pUndefMacroShadowItem a pL = 
            let l_arrayName = aName a
            in  "#undef " ++ pShowArrayTerm l_arrayName pL ++ breakline

{-
pShowKernel :: String -> PKernel -> String
pShowKernel l_name l_kernel = "Pochoir_Kernel_" ++ show dim ++ "D(" ++ l_name ++ ", " ++
    pShowKernelParams (kParams l_kernel) ++ ")" ++ show (kStmt l_kernel) ++
    breakline ++ "Pochoir_Kernel_end" ++ breakline
        where dim = length (kParams l_kernel) - 1
-}
pShowKernelParams :: [String] -> String
pShowKernelParams l_kernel_params = intercalate ", " l_kernel_params

pShowArrayGaps :: Int -> [PArray] -> String
pShowArrayGaps _ [] = ""
pShowArrayGaps l_rank l_array = breakline ++ "int " ++ 
        intercalate ", " (map (getArrayGaps (l_rank-1)) l_array) ++ ";"

{-
pShowInteriorKernel :: String -> PKernel -> String
pShowInteriorKernel l_name l_kernel =
    let l_iter = kIter l_kernel
        l_interiorStmts = transStmts (kStmt l_kernel)
                            $ transInterior $ getArrayName . unionArrayIter $ l_iter 
        l_interiorKernel = PKernel { kName = l_name,
                                     kParams = kParams l_kernel,
                                     kStmt = l_interiorStmts,
                                     kIter = kIter l_kernel 
                                   }
    in  pShowAutoKernel l_name l_interiorKernel

pShowTypeKernel :: [PArray] -> String -> PKernel -> String
pShowTypeKernel l_sArrayInUse l_name l_kernel =
    let l_iter = kIter l_kernel
    --  in  pShowShadowArrayInUse l_array ++ pShowAutoKernel l_name l_kernel 
    in  pShowShadowArrayInUse l_sArrayInUse ++ pShowAutoKernel l_name l_kernel 
-}

-- AutoKernel is a de-sugared kernel
pShowAutoKernel :: String -> PKernel -> String
pShowAutoKernel l_name l_kernel = 
    let l_params = zipWith (++) (repeat "int ") (kParams l_kernel)
    in  "/* known! */ auto " ++ l_name ++ " = [&] (" ++ 
        pShowKernelParams l_params ++ ") {" ++ breakline ++
        show (kStmt l_kernel) ++
        breakline ++ "};" ++ breakline

pShowUnrolledMacroKernels :: String -> String -> PStencil -> [PName] -> String
pShowUnrolledMacroKernels _ _ _ [] = ""
pShowUnrolledMacroKernels l_macro l_name l_stencil l_kL@(l_kernel:l_kernels) = 
    let l_t = "t"
        l_rank = sRank l_stencil
    in  breakline ++ "/* known! */ auto " ++ l_name ++ " = [&] (" ++
        "int t0, int t1, grid_info<" ++ show l_rank ++ "> const & grid) {" ++
        breakline ++ "grid_info<" ++ show l_rank ++ "> l_grid = grid;" ++
        breakline ++ "grid_info<" ++ show l_rank ++ "> l_phys_grid = " ++ 
        sName l_stencil ++ ".get_phys_grid();" ++
        breakline ++ pShowTimeLoopHeader l_t ++ breakline ++
        pShowSingleKernelIteration l_t l_macro l_rank l_kL ++
        breakline ++ pShowTimeLoopTail ++ breakline ++ "};\n"
    
pShowSingleKernelIteration :: String -> String -> Int -> [PName] -> String
pShowSingleKernelIteration _ _ _ [] = ""
pShowSingleKernelIteration l_t l_macro l_rank (l_kernel:l_kernels) =
    breakline ++ "meta_grid_" ++ l_macro ++ "<" ++ show l_rank ++ ">::single_step(" ++
    l_t ++ ", l_grid, l_phys_grid, " ++ l_kernel ++ ");" ++ 
    breakline ++ pAdjustTrape l_rank ++ breakline ++ pAdjustT ++ 
    breakline ++ pShowSingleKernelIteration l_t l_macro l_rank l_kernels
        where pAdjustT = if null l_kernels then "" else "++t;" 

pShowMacroKernel :: PName -> PKernel -> String
pShowMacroKernel l_macro l_kernel =
    let l_iters = kIter l_kernel
        l_name = l_macro ++ "_" ++ kName l_kernel
        l_sArrayInUse = unionArrayIter l_iters
        shadowArrayInUse = pDefMacroArrayInUse l_macro l_sArrayInUse (kParams l_kernel)
        unshadowArrayInUse = pUndefMacroArrayInUse l_sArrayInUse (kParams l_kernel)
    in  shadowArrayInUse ++ pShowAutoKernel l_name l_kernel ++ unshadowArrayInUse

pShowObaseKernel :: String -> PKernel -> String
pShowObaseKernel l_name l_kernel = 
    let l_rank = length (kParams l_kernel) - 1
        l_iter = kIter l_kernel
        l_array = unionArrayIter l_iter
        l_t = head $ kParams l_kernel
    in  breakline ++ "auto " ++ l_name ++ " = [&] (" ++
        "int t0, int t1, grid_info<" ++ show l_rank ++ "> const & grid) {" ++ 
        breakline ++ "grid_info<" ++ show l_rank ++ "> l_grid = grid;" ++
        pShowArrayGaps l_rank l_array ++
        breakline ++ pShowRankAttr l_rank "stride" l_array ++ breakline ++
        pShowTimeLoopHeader l_t ++ breakline ++
        pShowIterSet l_iter (kParams l_kernel)++
        breakline ++ pShowObaseForHeader l_rank l_iter (tail $ kParams l_kernel) ++
        breakline ++ pShowObaseStmt l_kernel ++ breakline ++ pShowObaseForTail l_rank ++
        pAdjustTrape l_rank ++ breakline ++ pShowTimeLoopTail ++ breakline ++ "};\n"

pShowPointerKernel :: String -> PKernel -> String
pShowPointerKernel l_name l_kernel = 
    let l_rank = length (kParams l_kernel) - 1
        l_iter = kIter l_kernel
        l_array = unionArrayIter l_iter
        l_t = head $ kParams l_kernel
    in  breakline ++ "auto " ++ l_name ++ " = [&] (" ++
        "int t0, int t1, grid_info<" ++ show l_rank ++ "> const & grid) {" ++ 
        breakline ++ "grid_info<" ++ show l_rank ++ "> l_grid = grid;" ++
        pShowPointers l_iter ++ breakline ++ 
        pShowArrayInfo l_array ++ pShowArrayGaps l_rank l_array ++
        breakline ++ pShowRankAttr l_rank "stride" l_array ++ breakline ++
        pShowTimeLoopHeader l_t ++ breakline ++
        pShowPointerSet l_iter (kParams l_kernel)++
        breakline ++ pShowPointerForHeader l_rank True l_iter (tail $ kParams l_kernel) ++
        breakline ++ pShowPointerStmt True l_kernel ++ breakline ++ pShowObaseForTail l_rank ++
        pAdjustTrape l_rank ++ breakline ++ pShowTimeLoopTail ++ breakline ++ "};\n"

pShowCachingKernel :: String -> PKernel -> String
pShowCachingKernel l_name l_kernel = 
    let l_rank = length (kParams l_kernel) - 1
        l_iters = kIter l_kernel
        l_stmts = kStmt l_kernel
        l_rdIters = pGetReadIters l_iters
        l_wrIters = pGetWriteIters l_iters 
        l_arrays = unionArrayIter l_iters
        l_t = head $ kParams l_kernel
    in  breakline ++ "auto " ++ l_name ++ " = [&] (" ++
        "int t0, int t1, grid_info<" ++ show l_rank ++ "> const & grid) {" ++ 
        breakline ++ "grid_info<" ++ show l_rank ++ "> l_grid = grid;" ++
        pShowPointers l_iters ++ breakline ++ 
        pShowArrayInfo l_arrays ++ pShowArrayGaps l_rank l_arrays ++
        breakline ++ pShowRankAttr l_rank "stride" l_arrays ++ breakline ++
        pShowDeltaT ++ breakline ++ pAllocStack l_arrays ++ breakline ++ 
        pShowLocalCacheAttr l_rank l_arrays ++ breakline ++ 
        pShowLocalCopyIn l_rank (kParams l_kernel) l_arrays l_rdIters ++ breakline ++ 
        pShowLocalGridForComp l_rank (head l_arrays) ++ breakline ++
        pShowTimeLoopHeader l_t ++ breakline ++
        pShowPointerSetLocal l_iters (kParams l_kernel) ++
        breakline ++ pShowPointerForHeader l_rank False l_iters (tail $ kParams l_kernel) ++
        breakline ++ pShowPointerStmt False l_kernel ++ breakline ++ 
        pShowObaseForTail l_rank ++ pAdjustTrape l_rank ++ breakline ++ 
        pShowTimeLoopTail ++ breakline ++ 
        pShowLocalCopyOut l_rank (kParams l_kernel) l_arrays l_wrIters ++ breakline ++ 
        pFreeStack l_arrays ++ breakline ++ 
        "};\n"

-- we assume all Pochoir_Array associated with the same Pochoir object has the same size,
-- otherwise we will trigger an error at run-time!
pShowLocalGridForCopyOut :: Int -> String
pShowLocalGridForCopyOut l_rank =
    let l_offsets = pGenRankListEmpty l_rank 
        l_x0s = pShowXs 0 l_rank "" l_offsets ++ breakline
        l_x1s = pShowXs 1 l_rank "" l_offsets ++ breakline
    in  l_x0s ++ l_x1s

pGenRankListEmpty :: Int -> [String]
pGenRankListEmpty 1 = [""]
pGenRankListEmpty l_rank = "":pGenRankListEmpty (l_rank-1)

pShowLocalGridForCopyIn :: Int -> String
pShowLocalGridForCopyIn l_rank =
    let l_offsets = pGenRankList "l_dx_" l_rank
        l_x0s = pShowXs 0 l_rank " - " l_offsets ++ breakline
        l_x1s = pShowXs 1 l_rank " + " l_offsets ++ breakline
    in  l_x0s ++ l_x1s

pShowLocalGridForComp :: Int -> PArray -> String
pShowLocalGridForComp l_rank l_array =
    let l_a = aName l_array
        l_offsets = pGenRankList ("lc_" ++ l_a ++ "_begin_") l_rank
        l_x0s = pShowXs 0 l_rank " - " l_offsets ++ breakline
        l_x1s = pShowXs 1 l_rank " - " l_offsets ++ breakline
    in  l_x0s ++ l_x1s

pShowXs :: Int -> Int -> String -> [String] -> String
pShowXs l_end l_rank l_op l_offsets =
    let lc_xs = pGenRankListFull ("l_grid.x" ++ show l_end ++ "[") l_rank "]"
        l_xs = pGenRankListFull ("grid.x" ++ show l_end ++ "[") l_rank "]"
        l_rs = zipWith (pIns l_op) l_xs l_offsets
        l_rrs = map (flip (++) ";") l_rs
        -- l_r1s = zipWith (pIns " + ") l_r0s l_dx0s
    in  intercalate breakline $ zipWith (pIns " = ") lc_xs l_rrs

-- PName : list of kernel parameters
pShowPointerSetLocal :: [Iter] -> [PName] -> String
pShowPointerSetLocal [] _ = ""
pShowPointerSetLocal iL@(i:is) l_kernelParams = concatMap pShowPointerSetLocalTerm iL
    where pShowPointerSetLocalTerm (iterName, array, dim, rw) = 
            let l_arrayName = aName array
                l_rank = length l_kernelParams - 1
                l_arrayBaseName = "lc_" ++ l_arrayName
                l_arrayTotalSize = "lc_" ++ l_arrayName ++ "_total_size"
                l_arrayStrideList = pGenRankList ("lc_" ++ l_arrayName ++ "_stride_") l_rank 
                l_transDimList = tail $ pShowTransDim dim l_kernelParams
                -- l_begins = pGenRankList ("lc_" ++ l_arrayName ++ "_begin_") l_rank
                -- l_offsets = zipWith (pInsParens " - ") (map show l_transDimList) l_begins
                l_arraySpaceOffset = 
                    intercalate " + " $ zipWith (pIns " * ") (map show l_transDimList) l_arrayStrideList
                l_arrayTimeOffset = (pShowTimeOffset (aToggle array) (head dim)) ++ 
                                    " * " ++ l_arrayTotalSize
            in  breakline ++ iterName ++ " = " ++ l_arrayBaseName ++ " + " ++ 
                l_arrayTimeOffset ++ " + " ++ l_arraySpaceOffset ++ ";" 

pShowLocalCopyOut :: Int -> [PName] -> [PArray] -> [Iter] -> String
pShowLocalCopyOut l_rank l_kParams l_arrays l_wrIters =
    "/* Copy Out */" ++ breakline ++ 
    pShowLocalGridForCopyOut l_rank ++ breakline ++
    pShowTimeLoopHeader (head l_kParams) ++ breakline ++
    (intercalate breakline $ map (pShowCopyInOutBase False l_rank $ head l_kParams) l_wrIters) ++ breakline ++
    pShowCopyLoopHeader False l_rank (tail l_kParams) l_arrays ++ breakline ++
    (intercalate breakline $ map (pShowCopyInOutBody False) l_arrays) ++ breakline ++
    pShowCopyLoopTail l_rank ++ breakline ++ pAdjustTrape l_rank ++ breakline ++ 
    pShowTimeLoopTail ++ breakline 

pShowLocalCopyIn :: Int -> [PName] -> [PArray] -> [Iter] -> String
pShowLocalCopyIn l_rank l_kParams l_arrays l_rdIters =
    "/* Copy In */" ++ breakline ++ 
    pShowLocalGridForCopyIn l_rank ++ breakline ++
    pShowTimeLoopHeader (head l_kParams) ++ breakline ++
    (intercalate breakline $ map (pShowCopyInOutBase True l_rank $ head l_kParams) l_rdIters) ++ breakline ++
    pShowCopyLoopHeader True l_rank (tail l_kParams) l_arrays ++ breakline ++
    (intercalate breakline $ map (pShowCopyInOutBody True) l_arrays) ++ breakline ++
    pShowCopyLoopTail l_rank ++ breakline ++ pAdjustTrape l_rank ++ breakline ++
    pShowTimeLoopTail ++ breakline

pShowCopyInOutBody :: Bool -> PArray -> String
pShowCopyInOutBody l_inOut l_array =
    let l_a = aName l_array
        l_suffix = if l_inOut == True then "_in" else "_out"
        lc_iter = "lc_" ++ l_a ++ l_suffix ++ "[0]"
        l_iter = "l_" ++ l_a ++ l_suffix ++ "[0]"
    in  if l_inOut == True 
            then lc_iter ++ " = " ++ l_iter ++ ";"
            else l_iter ++ " = " ++ lc_iter ++ ";"

pShowCopyLoopHeader :: Bool -> Int -> [PName] -> [PArray] -> String
pShowCopyLoopHeader l_inOut 1 l_kSpatialParams l_arrays =
    let l_loopVar = head l_kSpatialParams
        lc_incGap = intercalate ", " $ map (pShowIncGap l_inOut "lc_" 1) l_arrays
        l_incGap = intercalate ", " $ map (pShowIncGap l_inOut "l_" 1) l_arrays
    in  "for (int " ++ 
        l_loopVar ++ " = l_grid.x0[" ++ show (1-1) ++ "]; " ++ 
        l_loopVar ++ " < l_grid.x1[" ++ show (1-1) ++ "]; ++" ++
        l_loopVar ++ ", " ++ breakline ++ lc_incGap ++ ", " ++ 
        breakline ++ l_incGap ++ ") {"
pShowCopyLoopHeader l_inOut l_rank l_kSpatialParams l_arrays =
    let lc_gaps = intercalate breakline $ map (pShowArrayGap "lc_" l_rank) l_arrays
        l_gaps = intercalate breakline $ map (pShowArrayGap "l_" l_rank) l_arrays
        l_loopVar = head l_kSpatialParams
        lc_incGap = intercalate ", " $ map (pShowIncGap l_inOut "lc_" l_rank) l_arrays
        l_incGap = intercalate ", " $ map (pShowIncGap l_inOut "l_" l_rank) l_arrays
    in  lc_gaps ++ breakline ++ l_gaps ++ breakline ++ "for (int " ++ 
        l_loopVar ++ " = l_grid.x0[" ++ show (l_rank-1) ++ "]; " ++ 
        l_loopVar ++ " < l_grid.x1[" ++ show (l_rank-1) ++ "]; ++" ++
        l_loopVar ++ ", " ++ breakline ++ lc_incGap ++ ", " ++ 
        breakline ++ l_incGap ++ ") {" ++
        breakline ++ pShowCopyLoopHeader l_inOut (l_rank-1) (tail l_kSpatialParams) l_arrays

pShowIncGap :: Bool -> String -> Int -> PArray -> String
pShowIncGap l_inOut pre 1 l_array =
    let l_a = aName l_array
        l_suffix = if l_inOut == True then "_in" else "_out"
        l_iter = pre ++ l_a ++ l_suffix
    in  "++" ++ l_iter
pShowIncGap l_inOut pre l_rank l_array =
    let l_a = aName l_array
        l_suffix = if l_inOut == True then "_in" else "_out"
        l_iter = pre ++ l_a ++ l_suffix
        l_gap = pre ++ "gap_" ++ l_a ++ "_" ++ show (l_rank-1)
    in  l_iter ++ " += " ++ l_gap

pShowArrayGap :: String -> Int -> PArray -> String
pShowArrayGap pre l_rank l_array =
    let l_a = aName l_array
        l_gap = pre ++ "gap_" ++ l_a ++ "_" ++ show (l_rank-1)
        l_stride = pre ++ l_a ++ "_stride_" ++ show (l_rank-1)
        l_pre_stride = pre ++ l_a ++ "_stride_" ++ show (l_rank-2)
        l_offset = "(l_grid.x0[" ++ show (l_rank-2) ++ "] - l_grid.x1[" ++ show (l_rank-2) ++ "])"
    in  "int " ++ l_gap ++ " = " ++ l_stride ++ " + " ++ l_offset ++ " * " ++ l_pre_stride ++ ";"

pShowCopyLoopTail :: Int -> String
pShowCopyLoopTail 1 = "}"
pShowCopyLoopTail l_rank = "}" ++ breakline ++ pShowCopyLoopTail (l_rank-1)

pShowCopyInOutBase :: Bool -> Int -> PName -> Iter -> String
pShowCopyInOutBase l_inOut l_rank l_t l_rwIter =
    let l_array = pIterArray l_rwIter 
        l_type = show (aType l_array) ++ " * "
        l_a = aName l_array
        lc_base = "lc_" ++ l_a
        l_base = l_a ++ "_base"
        l_suffix = if l_inOut == True then "_in" else "_out"
        l_t_dim = head $ pIterDims l_rwIter 
        lc_iter = "lc_" ++ l_a ++ l_suffix
        l_iter = "l_" ++ l_a ++ l_suffix
        l_toggle = aToggle l_array
        lc_total_size = "lc_" ++ l_a ++ "_total_size"
        l_total_size = "l_" ++ l_a ++ "_total_size"
        lc_strides = pGenRankList ("lc_" ++ l_a ++ "_stride_") l_rank
        l_strides = pGenRankList ("l_" ++ l_a ++ "_stride_") l_rank
        l_grid_begins = pGenRankListFull "l_grid.x0[" l_rank "]"
        lc_begins = pGenRankList ("lc_" ++ l_a ++ "_begin_") l_rank
        lc_grid_offsets = zipWith (pInsParens " - ") l_grid_begins lc_begins
        lc_offsets = zipWith (pIns " * ") lc_grid_offsets lc_strides
        lc_a_base = l_type ++ lc_iter ++ " = " ++ lc_base ++ " + " ++
                  (pShowTimeOffset l_toggle l_t_dim) ++
                  " * " ++ lc_total_size ++ " + " ++ 
                  intercalate " + " lc_offsets ++ ";" ++ breakline
        l_grid_offsets = zipWith (pIns " * ") l_grid_begins l_strides
        l_a_base = l_type ++ l_iter ++ " = " ++ l_base ++ " + " ++
                 (pShowTimeOffset l_toggle l_t_dim) ++
                 " * " ++ l_total_size ++ " + " ++
                 intercalate " + " l_grid_offsets ++ ";" ++ breakline
    in  lc_a_base ++ l_a_base

pShowTimeLoopHeader :: String -> String
pShowTimeLoopHeader l_t =
    "for (int " ++ l_t ++ " = t0; " ++ l_t ++ " < t1; ++" ++ l_t ++ ") {"

pShowTimeLoopTail :: String
pShowTimeLoopTail = "} /* end for t */"

pShowCopyInLoop :: Int -> [PName] -> PArray -> String
pShowCopyInLoop l_rank l_kSpatialParams l_array =
    let l_forHead = pShowLocalForHeader l_rank l_kSpatialParams l_array
        l_forTail = pShowLocalForTail l_rank
        l_toggle = aToggle l_array
        l_body = pShowLocalCopyInBody l_toggle l_kSpatialParams l_array
    in  l_forHead ++ l_body ++ l_forTail

pShowLocalCopyInBody :: Int -> [PName] -> PArray -> String
pShowLocalCopyInBody 1 l_kSpatialParams l_array =
    let l_a = aName l_array
        l_rank = aRank l_array
        l_local_addr = "lc_" ++ l_a ++ "_in_" ++ show 0
        l_global_addr = "l_" ++ l_a ++ "_in_" ++ show 0
        l_local_offset = intercalate " + " $ zipWith (pIns " * ") l_kSpatialParams (pGenRankList ("lc_" ++ l_a ++ "_stride_") l_rank)
        l_global_offset = intercalate " + " $ zipWith (pIns " * ") l_kSpatialParams (pGenRankList ("l_" ++ l_a ++ "_stride_") l_rank)
    in  l_local_addr ++ "[" ++ l_local_offset ++ "] = " ++ l_global_addr ++ "[" ++ l_global_offset ++ "];" ++ breakline
pShowLocalCopyInBody l_toggle l_kSpatialParams l_array =
    let l_a = aName l_array
        l_rank = aRank l_array
        l_local_addr = "lc_" ++ l_a ++ "_in_" ++ show (l_toggle-1)
        l_global_addr = "l_" ++ l_a ++ "_in_" ++ show (l_toggle-1)
        l_local_offset = intercalate " + " $ zipWith (pIns " * ") l_kSpatialParams (pGenRankList ("lc_" ++ l_a ++ "_stride_") l_rank)
        l_global_offset = intercalate " + " $ zipWith (pIns " * ") l_kSpatialParams (pGenRankList ("l_" ++ l_a ++ "_stride_") l_rank)
    in  l_local_addr ++ "[" ++ l_local_offset ++ "] = " ++ l_global_addr ++ "[" ++ l_global_offset ++ "];" ++ breakline ++ pShowLocalCopyInBody (l_toggle-1) l_kSpatialParams l_array

pIns :: String -> String -> String -> String
pIns mid a b = a ++ mid ++ b

pInsParens :: String -> String -> String -> String
pInsParens mid a b = "(" ++ a ++ mid ++ b ++ ")"

pGenRankList :: String -> Int -> [String]
pGenRankList pre 1 = [pre ++ show 0]
pGenRankList pre l_rank = [pre ++ show (l_rank-1)] ++ pGenRankList pre (l_rank-1)

pGenRankListFull :: String -> Int -> String -> [String]
pGenRankListFull pre 1 suf = [pre ++ show 0 ++ suf]
pGenRankListFull pre l_rank suf = [pre ++ show (l_rank-1) ++ suf] ++ 
                                  pGenRankListFull pre (l_rank-1) suf

pShowLocalForHeader :: Int -> [PName] -> PArray -> String
pShowLocalForHeader 1 l_kSpatialParams l_array =
    let l_loopVar = head l_kSpatialParams
        l_begin = show 0
        l_end = "lc_" ++ aName l_array ++ "_size_" ++ show 0
    in  "for (int " ++ l_loopVar ++ " = " ++ l_begin ++ "; " ++ 
        l_loopVar ++ " < " ++ l_end ++ "; " ++ "++" ++ l_loopVar ++ ") {" ++ 
        breakline
pShowLocalForHeader n l_kSpatialParams l_array =
    let l_loopVar = head l_kSpatialParams
        l_begin = show 0 
        l_end = "lc_" ++ aName l_array ++ "_size_" ++ show (n-1)
    in  "for (int " ++ l_loopVar ++ " = " ++ l_begin ++ "; " ++ 
        l_loopVar ++ " < " ++ l_end ++ "; " ++ "++" ++ l_loopVar ++ ") {" ++ 
        breakline ++ pShowLocalForHeader (n-1) (tail l_kSpatialParams) l_array

pShowLocalForTail :: Int -> String
pShowLocalForTail 1 = "}"
pShowLocalForTail n = "}" ++ pShowLocalForTail (n-1) 

pShowDeltaT :: String
pShowDeltaT = "const int lt = t1 - t0;"

pFreeHeap :: [PArray] -> String
pFreeHeap l_arrays = concatMap pFreeHeapItem l_arrays
    where pFreeHeapItem l_array =
            let l_a = aName l_array
                l_cache = "lc_" ++ l_a
            in  "free(" ++ l_cache ++ ");" ++ breakline

pAllocHeap :: [PArray] -> String
pAllocHeap l_arrays = concatMap pAllocHeapItem l_arrays
    where pAllocHeapItem l_array = 
            let l_type = show (aType l_array) 
                l_a = aName l_array
                l_size = " 2 * 120 * 120 "
                l_cache = " lc_" ++ l_a
            in  l_type ++ " *" ++ l_cache ++ " = new " ++ 
                l_type ++ " [" ++ l_size ++ "];" ++ breakline

pFreeStack :: [PArray] -> String
pFreeStack l_arrays = ""

pAllocStack :: [PArray] -> String
pAllocStack l_arrays = concatMap pAllocStackItem l_arrays
    where pAllocStackItem l_array = 
            let l_type = show (aType l_array) 
                l_a = aName l_array
                l_toggle = aToggle l_array
                l_size = show l_toggle ++ " * 120 * 120"
                l_cache = " lc_" ++ l_a
            in  l_type ++ l_cache ++ " [ " ++ l_size ++ " ];" ++ breakline

pShowLocalCacheAttr :: Int -> [PArray] -> String
pShowLocalCacheAttr l_rank l_arrays = 
    pShowColor l_rank ++ breakline ++ 
    pShowSlopes l_rank ++ breakline ++ 
    (intercalate breakline $ map (pShowEndIndex l_rank 0) l_arrays) ++ breakline ++ 
    (intercalate breakline $ map (pShowEndIndex l_rank 1) l_arrays) ++ breakline ++
    (intercalate breakline $ map (pShowLocalArraySize l_rank) l_arrays) ++ breakline ++
    (intercalate breakline $ map (pShowLocalStride l_rank) l_arrays) ++ breakline ++
    (intercalate breakline $ map (pShowLocalTotalSize l_rank) l_arrays) ++ breakline 

pShowLocalTotalSize :: Int -> PArray -> String
pShowLocalTotalSize  l_rank l_array =
    let l_a = aName l_array
        l_var = "lc_" ++ l_a ++ "_total_size" 
        l_total_size = pShowSizeSum "lc" l_rank l_array
    in  "const int " ++ l_var ++ " = " ++ l_total_size ++ ";"

pShowLocalStride :: Int -> PArray -> String
pShowLocalStride l_rank l_array = 
    "const int " ++ pShowLocalStrideTerm l_rank l_array

pShowLocalStrideTerm :: Int -> PArray -> String
pShowLocalStrideTerm 1 l_array = 
    let l_a = aName l_array
        l_dim = show 0
        l_strideName = "lc_" ++ l_a ++ "_stride_" ++ l_dim
        l_size = show 1
    in  l_strideName ++ " = " ++ l_size ++ ";"
pShowLocalStrideTerm n l_array =
    let l_a = aName l_array
        l_dim = show (n-1)
        l_strideName = "lc_" ++ l_a ++ "_stride_" ++ l_dim
        l_size = pShowSizeSum "lc" (n-1) l_array
    in  l_strideName ++ " = " ++ l_size ++ ", " ++ pShowLocalStrideTerm (n-1) l_array

pShowSizeSum :: String -> Int -> PArray -> String
pShowSizeSum pre 1 l_array =
    let l_a = aName l_array
        l_dim = show 0
        l_size = pre ++ "_" ++ l_a ++ "_size_" ++ l_dim
    in  l_size 
pShowSizeSum pre n l_array =
    let l_a = aName l_array
        l_dim = show (n-1) 
        l_size = pre ++ "_" ++ l_a ++ "_size_" ++ l_dim
    in  l_size ++ " * " ++ pShowSizeSum pre (n-1) l_array 

pShowLocalArraySize :: Int -> PArray -> String
pShowLocalArraySize 1 l_array =
    let l_a = aName l_array
        l_dim = show 0
        l_var = "lc_" ++ l_a ++ "_size_" ++ l_dim
        l_end = "lc_" ++ l_a ++ "_end_" ++ l_dim
        l_begin = "lc_" ++ l_a ++ "_begin_" ++ l_dim
    in  "const int " ++ l_var ++ " = " ++ l_end ++ " - " ++ l_begin ++ " + 1;"
pShowLocalArraySize n l_array =
    let l_a = aName l_array
        l_dim = show (n-1) 
        l_var = "lc_" ++ l_a ++ "_size_" ++ l_dim
        l_end = "lc_" ++ l_a ++ "_end_" ++ l_dim
        l_begin = "lc_" ++ l_a ++ "_begin_" ++ l_dim
    in  "const int " ++ l_var ++ " = " ++ l_end ++ " - " ++ l_begin ++ " + 1;" ++ 
        breakline ++ pShowLocalArraySize (n-1) l_array

pShowSlopes :: Int -> String
pShowSlopes l_rank = 
    let l_dxs = pGenRankList "const int l_dx_" l_rank
        l_slopes = pGenRankListFull "a.slope(" l_rank ")"
        l_rs = map (flip (++) ";")  l_slopes
    in  intercalate breakline $ zipWith (pIns " = ") l_dxs l_rs

pShowColor :: Int -> String
pShowColor 1 = "const bool black_0 = (grid.dx0[0] >= 0 & grid.dx1[0] <= 0);"
pShowColor n = 
    let l_dim = show (n - 1)
        l_begin_slope = "grid.dx0[" ++ l_dim ++ "]"
        l_end_slope = "grid.dx1[" ++ l_dim ++ "]"
        l_var = "black_" ++ l_dim
    in  "const bool " ++ l_var ++ " = (" ++ l_begin_slope ++ " >= 0 & " ++
        l_end_slope ++ " <= 0);" ++ breakline ++ pShowColor (n-1)

pShowEndIndex :: Int -> Int -> PArray -> String
pShowEndIndex l_rank l_end l_array =
    let l_a = aName l_array
        l_endVar = if l_end == 0 then "_begin_" else "_end_"
        l_type = "const int "
        l_lvalues = pGenRankList ("lc_" ++ l_a ++ l_endVar) l_rank
        l_lefts = map ((++) l_type) l_lvalues
        l_blacks = pGenRankList "black_" l_rank
        l_slopes = pGenRankList "l_dx_" l_rank
        l_grid_xs = if l_end == 0 
                        then pGenRankListFull "grid.x0[" l_rank "]"
                        else pGenRankListFull "grid.x1[" l_rank "]"
        l_grid_dxs = if l_end == 0
                        then pGenRankListFull "grid.dx0[" l_rank "]"
                        else pGenRankListFull "grid.dx1[" l_rank "]"
        l_glue = if l_end == 0 then " - " else " + "
        l_r1s = zipWith (pIns l_glue) l_grid_xs l_slopes
        l_r20s = map (flip (++) " * lt") l_grid_dxs
        l_r21s = zipWith (pIns l_glue) l_r20s l_slopes
        l_r2s = zipWith (pIns " + ") l_grid_xs l_r21s
        l_rr0s = zipWith (pIns " ? ") l_blacks l_r1s
        l_rr1s = zipWith (pIns " : ") l_rr0s l_r2s
        l_rights = map (flip (++) ";") l_rr1s
    in  intercalate breakline $ zipWith (pIns " = ") l_lefts l_rights

pShowOptPointerKernel :: String -> PKernel -> String
pShowOptPointerKernel l_name l_kernel = 
    let l_rank = length (kParams l_kernel) - 1
        l_iter = kIter l_kernel
        l_array = unionArrayIter l_iter
        l_t = head $ kParams l_kernel
    in  breakline ++ "auto " ++ l_name ++ " = [&] (" ++
        "int t0, int t1, grid_info<" ++ show l_rank ++ "> const & grid) {" ++ 
        breakline ++ "grid_info<" ++ show l_rank ++ "> l_grid = grid;" ++
        pShowPointers l_iter ++ breakline ++ 
        pShowArrayInfo l_array ++ pShowArrayGaps l_rank l_array ++
        breakline ++ pShowRankAttr l_rank "stride" l_array ++ breakline ++
        pShowTimeLoopHeader l_t ++ breakline ++
        pShowOptPointerSet l_iter (kParams l_kernel)++
        breakline ++ pShowPointerForHeader l_rank True l_iter (tail $ kParams l_kernel) ++
        breakline ++ pShowOptPointerStmt l_kernel ++ breakline ++ pShowObaseForTail l_rank ++
        pAdjustTrape l_rank ++ breakline ++ pShowTimeLoopTail ++ breakline ++ "};\n"

pShowCPointerKernel :: String -> PKernel -> String
pShowCPointerKernel l_name l_kernel = 
    let l_rank = length (kParams l_kernel) - 1
        l_iter = kIter l_kernel
        l_array = unionArrayIter l_iter
        l_t = head $ kParams l_kernel
    in  breakline ++ "auto " ++ l_name ++ " = [&] (" ++
        "int t0, int t1, grid_info<" ++ show l_rank ++ "> const & grid) {" ++ 
        breakline ++ "grid_info<" ++ show l_rank ++ "> l_grid = grid;" ++
        pShowArrayInfo l_array ++ 
        breakline ++ pShowRankAttr l_rank "stride" l_array ++ breakline ++
        pShowRefMacro (kParams l_kernel) l_array ++
        "for (int " ++ l_t ++ " = t0; " ++ l_t ++ " < t1; ++" ++ l_t ++ ") { " ++ 
        pShowTimeLoopHeader l_t ++ breakline ++
        breakline ++ pShowRawForHeader (tail $ kParams l_kernel) ++
        breakline ++ pShowCPointerStmt l_kernel ++ breakline ++ pShowObaseForTail l_rank ++
        pAdjustTrape l_rank ++ breakline ++ pShowRefUnMacro l_array ++ pShowTimeLoopTail ++ 
        breakline ++ "};\n"

pShowCPointerStmt :: PKernel -> String
pShowCPointerStmt l_kernel = 
    let oldStmts = kStmt l_kernel
        l_iter = kIter l_kernel
        obaseStmts = transStmts oldStmts $ transCPointer l_iter
    in show obaseStmts

transCPointer :: [Iter] -> Expr -> Expr
transCPointer l_iters (PVAR q v dL) =
    case pIterLookup (v, dL) l_iters of
        Nothing -> PVAR q v dL
        Just iterName -> VAR q $ pRef v dL
transCPointer l_iters e = e

pRef :: PName -> [DimExpr] -> String
pRef a dL = "ref_" ++ a ++ "(" ++ (intercalate ", " $ map show dL) ++ ")"

pShowRawForHeader :: [PName] -> String
pShowRawForHeader [] = ""
pShowRawForHeader pL@(p:ps) = 
    let len_pL = length pL
        idx = p
        l_rank = len_pL-1
        l_pragma = if l_rank == 0 then pShowPragma else ""
    in  l_pragma ++ 
        breakline ++ "for (int " ++ idx ++ " = l_grid.x0[" ++ show l_rank ++
        "]; " ++ idx ++ " < l_grid.x1[" ++ show l_rank ++ "]; ++" ++ idx ++ ") {" ++
        pShowRawForHeader ps 
 
pShowPragma :: String
pShowPragma = breakline ++ "#pragma ivdep"
-- pShowPragma = "#pragma ivdep" ++ breakline ++ "#pragma simd"

pShowRefUnMacro :: [PArray] -> String
pShowRefUnMacro [] = ""
pShowRefUnMacro (a:as) = 
    let l_name = aName a
    in  "#undef ref_" ++ l_name ++ breakline ++ breakline

pShowRefMacro :: [PName] -> [PArray] -> String
pShowRefMacro _ [] = ""
pShowRefMacro l_kernelParams aL@(a:as) =
    let l_name = aName a
        l_t = head l_kernelParams
        l_dims = tail l_kernelParams
        l_toggle = aToggle a
        l_rank = aRank a
    in  "#define ref_" ++ l_name ++ "(" ++ pShowKernelParams l_kernelParams ++
        ") " ++ l_name ++ "_base[" ++ pShowTimeOffset l_toggle (DimVAR l_t) ++
        " * l_" ++ l_name ++ "_total_size + " ++ 
        (intercalate " + " $ zipWith pMul l_dims $ pGenRankList ("l_" ++ l_name ++ "_stride_") l_rank) ++ "]" ++
        breakline ++ breakline ++ pShowRefMacro l_kernelParams as

pMul :: String -> String -> String
pMul a b = "(" ++ a ++ ") * " ++ b

pShowArrayInfo :: [PArray] -> String
pShowArrayInfo [] = ""
pShowArrayInfo arrayInUse = foldr pShowArrayInfoItem "" arrayInUse
    where pShowArrayInfoItem l_arrayItem str =
            let l_type = aType l_arrayItem
                l_a = aName l_arrayItem
            in  str ++ breakline ++ show l_type ++ " * " ++ l_a ++ "_base"  ++ 
                " = " ++ l_a ++ ".data();" ++ breakline ++
                "const int " ++ "l_" ++ l_a ++ "_total_size = " ++ l_a ++
                ".total_size();" ++ breakline

pShowRankAttr :: Int -> String -> [PArray] -> String
pShowRankAttr n _ [] = ""
pShowRankAttr n attr aL@(a:as) = "const int " ++ getAttrs n aL ++ ";\n"
    where getAttrs n aL@(a:as) = intercalate ", " $ concatMap (getAttr n) aL
          getAttr 1 a = let r = 0 
                        in  ["l_" ++ (aName a) ++ "_" ++ attr ++ "_" ++ 
                             show r ++ " = " ++ (aName a) ++ "." ++ attr ++ 
                             "(" ++ show r ++ ")"]
          getAttr n a = let r = n-1
                        in  ["l_" ++ (aName a) ++ "_" ++ attr ++ "_" ++ 
                             show r ++ " = " ++ (aName a) ++ "." ++ attr ++
                             "(" ++ show r ++ ")"] ++ getAttr (n-1) a

pShowPointers :: [Iter] -> String
pShowPointers [] = ""
pShowPointers iL@(i:is) = foldr pShowPointer "" iL
    where pShowPointer (nameIter, arrayInUse, dL, rw) str =
                str ++ breakline ++ (show $ aType arrayInUse) ++ " * " ++ nameIter ++ ";"

pShowPointerStmt :: Bool -> PKernel -> String
pShowPointerStmt global l_kernel = 
    let oldStmts = kStmt l_kernel
        l_iter = kIter l_kernel
        obaseStmts = transStmts oldStmts $ transPointer global l_iter
    in show obaseStmts

pShowObaseStmt :: PKernel -> String
pShowObaseStmt l_kernel = 
    let oldStmts = kStmt l_kernel
        l_iter = kIter l_kernel
        obaseStmts = transStmts oldStmts $ transIter l_iter
    in show obaseStmts

pShowOptPointerStmt :: PKernel -> String
pShowOptPointerStmt l_kernel = 
    let oldStmts = kStmt l_kernel
        l_iter = kIter l_kernel
        obaseStmts = transStmts oldStmts $ transOptPointer l_iter
    in show obaseStmts

transOptPointer :: [Iter] -> Expr -> Expr
transOptPointer l_iters (PVAR q v dL) =
    case pIterLookup (v, dL) l_iters of
        Nothing -> PVAR q v dL
        Just iterName -> VAR q $ "(*" ++ iterName ++ ")"
transOptPointer l_iters e = e

transPointer :: Bool -> [Iter] -> Expr -> Expr
transPointer global l_iters (PVAR q v dL) =
    let pre = if global == True then "l_" else "lc_"
    in  case pPointerLookup (v, dL) l_iters of
            Nothing -> PVAR q v dL
            Just (iterName, arrayInUse, des, rw) -> 
                BVAR iterName de
                    where de = simplifyDimExpr naive_de
                          naive_de = foldr plusCombDimExpr x $ zipWith mulDimExpr strideL $ tail $ excludeDimExpr dL des 
                          strideL = pGenRankList (pre ++ aName arrayInUse ++ "_stride_") (aRank arrayInUse) 
                          x = (DimINT 0)
transPointer _ _ e = e

plusCombDimExpr :: DimExpr -> DimExpr -> DimExpr
plusCombDimExpr e1 e2 = DimDuo "+" e1 e2

mulDimExpr :: String -> DimExpr -> DimExpr
mulDimExpr stride dim = DimDuo "*" (DimVAR stride) dim

excludeDimExpr :: [DimExpr] -> [DimExpr] -> [DimExpr]
excludeDimExpr dL [] = dL
excludeDimExpr [] _ = []
excludeDimExpr (d:ds) (r:rs) = (excludeDimExprItem d r):(excludeDimExpr ds rs)
    where excludeDimExprItem (DimVAR v) r = if (DimVAR v) == r then DimINT 0 else (DimVAR v)
          excludeDimExprItem (DimDuo bop e1 e2) r 
            | (DimDuo bop e1 e2) == r = DimINT 0
            | e1 == r = DimParen (DimDuo bop (DimINT 0) e2)
            | e2 == r = DimParen (DimDuo bop e1 (DimINT 0))
            | otherwise = (DimDuo bop (excludeDimExprItem e1 r) (excludeDimExprItem e2 r))
          excludeDimExprItem (DimINT n) r = if (DimINT n) == r then DimINT 0 else (DimINT n)
          excludeDimExprItem (DimParen e) r 
            | (DimParen e) == r = DimINT 0
            | e == r = DimINT 0
            | otherwise = DimParen (excludeDimExprItem e r)

pPointerLookup :: (PName, [DimExpr]) -> [Iter] -> Maybe Iter
pPointerLookup (v, dL) [] = Nothing
pPointerLookup (v, dL) ((iterName, arrayInUse, dL', rw):is)
    | v == aName arrayInUse && head dL == head dL' = Just (iterName, arrayInUse, dL', rw)
    | otherwise = pPointerLookup (v, dL) is

transIter :: [Iter] -> Expr -> Expr
transIter l_iters (PVAR q v dL) =
    case pIterLookup (v, dL) l_iters of
        Nothing -> PVAR q v dL
        Just iterName -> VAR q iterName
transIter l_iters e = e

pShowIterSet :: [Iter] -> [PName] -> String
pShowIterSet iL@(i:is) l_kernelParams = concat $ map pShowIterSetTerm iL
    where pShowIterSetTerm (name, array, dim, rw) = 
            breakline ++ name ++ ".set(" ++ show (pShowTransDim dim l_kernelParams) ++ ");" 
            --
-- PName : list of kernel parameters
pShowOptPointerSet :: [Iter] -> [PName] -> String
pShowOptPointerSet [] _ = ""
pShowOptPointerSet iL@(i:is) l_kernelParams = 
    let baseIters = transIterN 0 $ getBaseIter l_kernelParams iL
    in  pShowPointers baseIters ++ (concat $ map pShowOptPointerSetTerm baseIters) ++ pShowNonBaseIters baseIters iL
        where pShowOptPointerSetTerm (iterName, array, dim, rw) = 
                let l_arrayName = aName array
                    l_arrayBaseName = l_arrayName ++ "_base"
                    l_arrayTotalSize = "l_" ++ l_arrayName ++ "_total_size"
                    l_arrayStrideList = 
                        pGenRankList ("l_" ++ l_arrayName ++ "_stride_") (length l_kernelParams - 1)
                    l_transDimList = tail $ pShowTransDim dim l_kernelParams
                    l_arraySpaceOffset = 
                        intercalate " + " $ zipWith pCombineDim l_transDimList l_arrayStrideList
                    l_arrayTimeOffset = (pShowTimeOffset (aToggle array) (head dim)) ++ 
                                        " * " ++ l_arrayTotalSize
                in  breakline ++ iterName ++ " = " ++ l_arrayBaseName ++ " + " ++ 
                    l_arrayTimeOffset ++ " + " ++ l_arraySpaceOffset ++ ";" 

pShowNonBaseIters :: [Iter] -> [Iter] -> String
pShowNonBaseIters _ [] = ""
pShowNonBaseIters bL iL@(i:is) = pShowShiftFromBase i bL ++ pShowNonBaseIters bL is
    where pShowShiftFromBase _ [] = "/* NO baseIter found! */"
          pShowShiftFromBase (iName, iArray, iDims, iRW) ((bName, bArray, bDims, bRW):bs)
            | iArray == bArray && head iDims == head bDims =
                let l_arrayName = aName iArray
                    l_arrayStrideList =
                        pGenRankList ("l_" ++ l_arrayName ++ "_stride_") (length iDims - 1) 
                    l_transDimList = map simplifyDimExpr $ zipWith excludeBaseDim (tail iDims) (tail bDims)
                    l_arraySpaceOffset = 
                        intercalate " + " $ zipWith pCombineDim l_transDimList l_arrayStrideList
                in  breakline ++ iName ++ " = "  ++ bName ++ " + " ++ l_arraySpaceOffset ++ ";"
            | otherwise = pShowShiftFromBase (iName, iArray, iDims, iRW) bs

excludeBaseDim :: DimExpr -> DimExpr -> DimExpr
excludeBaseDim (DimVAR i) b 
    | DimVAR i == b = (DimINT 0)
    | otherwise = DimVAR i
excludeBaseDim (DimINT n) b
    | DimINT n == b = (DimINT 0)
    | otherwise = DimINT n
excludeBaseDim (DimDuo bop e1 e2) b = 
    DimDuo bop (excludeBaseDim e1 b) (excludeBaseDim e2 b)
excludeBaseDim (DimParen e) b = DimParen (excludeBaseDim e b)

-- PName : list of kernel parameters
pShowPointerSet :: [Iter] -> [PName] -> String
pShowPointerSet [] _ = ""
pShowPointerSet iL@(i:is) l_kernelParams = concatMap pShowPointerSetTerm iL
    where pShowPointerSetTerm (iterName, array, dim, rw) = 
            let l_arrayName = aName array
                l_arrayBaseName = l_arrayName ++ "_base"
                l_arrayTotalSize = "l_" ++ l_arrayName ++ "_total_size"
                l_arrayStrideList = 
                    pGenRankList ("l_" ++ l_arrayName ++ "_stride_") (length l_kernelParams - 1) 
                l_transDimList = tail $ pShowTransDim dim l_kernelParams
                l_arraySpaceOffset = 
                    intercalate " + " $ zipWith pCombineDim l_transDimList l_arrayStrideList
                l_arrayTimeOffset = (pShowTimeOffset (aToggle array) (head dim)) ++ 
                                    " * " ++ l_arrayTotalSize
            in  breakline ++ iterName ++ " = " ++ l_arrayBaseName ++ " + " ++ 
                l_arrayTimeOffset ++ " + " ++ l_arraySpaceOffset ++ ";" 

pShowTimeOffset :: Int -> DimExpr -> String
pShowTimeOffset toggle tDim 
    | toggle == 2 = "((" ++ show tDim ++ ")" ++ " & 0x1" ++ ")"
    | toggle == 4 = "((" ++ show tDim ++ ")" ++ " & 0x11" ++ ")"
    | otherwise = "((" ++ show tDim ++ ") % " ++ show toggle ++ ")"

pCombineDim :: DimExpr -> String -> String
-- l_stride_pa_0 may NOT necessary be "1", 
-- plus that we have already set all strides to be of type "const int"
pCombineDim de stride = "(" ++ show de ++ ") * " ++ stride

pShowTransDim :: [DimExpr] -> [PName] -> [DimExpr]
pShowTransDim (d:ds) (p:ps) = 
    let l_rank = length ds
    in  (d:(pTransDim 1 l_rank ds ps))

pTransDim :: Int -> Int -> [DimExpr] -> [PName] -> [DimExpr]
pTransDim n r [] _ = []
pTransDim n r dL@(d:ds) pL@(p:ps) = pTransDimTerm n r d p : pTransDim (n+1) r ds ps
    where pTransDimTerm n r (DimVAR v) p
              | v == p = DimVAR ("l_grid.x0[" ++ show (r-n) ++ "]")
              | otherwise = DimVAR v
          pTransDimTerm n r (DimINT i) p = DimINT i
          pTransDimTerm n r (DimDuo bop e1 e2) p = 
              DimDuo bop (pTransDimTerm n r e1 p) (pTransDimTerm n r e2 p)
        
{-
pShowIters :: [Iter] -> String
pShowIters [] = ""
pShowIters ((l_name, l_array, l_dim, l_rw):is) = 
    let l_type = aType l_array
        l_rank = aRank l_array
        l_toggle = aToggle l_array
        l_arrayName = aName l_array
    in breakline ++ "Pochoir_Iterator<" ++ show l_type ++ ", " ++ show l_rank ++ 
       ", " ++ show l_toggle ++ "> " ++
       l_name ++ "(" ++ l_arrayName ++ ");" ++ pShowIters is   
-}

unionArrayIter :: [Iter] -> [PArray]
unionArrayIter [] = []
unionArrayIter iL@(i:is) = union (getArrayItem i) (unionArrayIter is)
    where getArrayItem (_, a, _, _) = [a]

getArrayIter :: [Iter] -> [PName]
getArrayIter [] = []
getArrayIter iL@(i:is) = (getArrayItem i) ++ (getArrayIter is)
    where getArrayItem (_, a, _, _) = [aName a]

pShowObaseForTail :: Int -> String
pShowObaseForTail n 
    | n == 0 = "/* end for (sub-trapezoid) */ "
    | otherwise = "} " ++ pShowObaseForTail (n-1)

pAdjustTrape :: Int -> String
pAdjustTrape l_rank = 
    breakline ++ "/* Adjust sub-trapezoid! */" ++
    breakline ++ "for (int i = 0; i < " ++ show l_rank ++ "; ++i) {" ++ 
    breakline ++ "\tl_grid.x0[i] += l_grid.dx0[i]; l_grid.x1[i] += l_grid.dx1[i];" ++
    breakline ++ "}" 

-- pL is the parameter list of original user supplied computing kernel
pShowObaseForHeader :: Int -> [Iter] -> [PName] -> String
pShowObaseForHeader _ _ [] = ""
pShowObaseForHeader 1 iL pL = 
                           breakline ++ pShowForHeader 0 True (unionArrayIter iL) pL ++ 
                           pShowIterComma iL ++
                           breakline ++ intercalate (", " ++ breakline) 
                                        (map ((++) "++" . getIterName) iL) ++ ") {"
pShowObaseForHeader n iL pL = 
                           breakline ++ pShowForHeader (n-1) True (unionArrayIter iL) pL ++ 
                           pShowIterComma iL ++
                           breakline ++ intercalate (", " ++ breakline) 
                                     (zipWith wrapIterInc
                                        (map (getArrayGap (n-1)) (getArrayIter iL))
                                        (map getIterName iL)) ++ 
                           ") {" ++ pShowObaseForHeader (n-1) iL pL
    where wrapIterInc gap iter = iter ++ ".inc(" ++ gap ++ ")"

-- pL is the parameter list of original user supplied computing kernel
pShowPointerForHeader :: Int -> Bool -> [Iter] -> [PName] -> String
pShowPointerForHeader _ _ _ [] = ""
pShowPointerForHeader 1 global iL pL = 
                           breakline ++ pShowPragma ++
                           breakline ++ pShowForHeader 0 global (unionArrayIter iL) pL ++  
                           pShowIterComma iL ++
                           breakline ++ intercalate (", " ++ breakline) 
                                        (map ((++) "++" . getIterName) iL) ++ ") {"
--                                        (map ((flip (++) "+=1") . getIterName) iL) ++ ") {"

pShowPointerForHeader n global iL pL = 
                           breakline ++ 
                           pShowForHeader (n-1) global (unionArrayIter iL) pL ++ 
                           pShowIterComma iL ++
                           breakline ++ intercalate (", " ++ breakline) 
                                     (zipWith wrapIterInc
                                        (map (getArrayGap (n-1)) (getArrayIter iL))
                                        (map getIterName iL)) ++ 
                           ") {" ++ pShowPointerForHeader (n-1) global iL pL
    where wrapIterInc gap iter = iter ++ " += " ++ gap 

pShowIterComma :: [Iter] -> String
pShowIterComma [] = ""
pShowIterComma iL@(i:is) = ", "

pShowForHeader :: Int -> Bool -> [PArray] -> [PName] -> String
pShowForHeader _ _ _ [] = ""
pShowForHeader i global aL pL = 
    let len_pL = length pL
        l_var = pL !! (len_pL - 1 - i)
        l_rank = show i
        l_begin = "l_grid.x0[" ++ l_rank ++ "]"
        l_end = "l_grid.x1[" ++ l_rank ++ "]"
    in  adjustGap i global aL ++ "for (int " ++ l_var ++ 
        " = " ++ l_begin ++ "; " ++ l_var ++ " < " ++ l_end ++ "; ++" ++ l_var
                    
adjustGap :: Int -> Bool -> [PArray] -> String
adjustGap l_rank _ [] = ""
adjustGap l_rank global aL@(a:as) = 
    if l_rank > 0 then pShowAdjustGap l_rank aL
             else ""
    where pShowAdjustGap l_rank [] = ""
          pShowAdjustGap l_rank aL@(a:as) = concatMap (pShowAdjustGapTerm l_rank) aL
          pShowAdjustGapTerm l_rank a = 
            let l_a = aName a
                l_gap = "gap_" ++ l_a ++ "_" ++ show l_rank
                pre = if global == True then "l_" else "lc_"
                l_stride = pre ++ l_a ++ "_stride_" ++ show l_rank 
                l_pre_stride = pre ++ l_a ++ "_stride_" ++ show (l_rank - 1)
                l_begin = "l_grid.x0[" ++ show (l_rank-1) ++ "]"
                l_end = "l_grid.x1[" ++ show (l_rank-1) ++ "]"
            in  l_gap ++ " = " ++ l_stride ++ " + (" ++ l_begin ++ " - " ++ l_end ++ 
                ") * " ++ l_pre_stride ++ ";" ++ breakline
          
getIterName :: Iter -> String
getIterName (name, _, _, _) = name

getArrayGaps :: Int -> PArray -> String
getArrayGaps 0 array = getArrayGap 0 (aName array)
getArrayGaps n array = getArrayGap n (aName array) ++ ", " ++ getArrayGaps (n-1) array

getArrayGap :: Int -> PName -> String
getArrayGap n array = "gap_" ++ array ++ "_" ++ show n

pShowShapes :: [[Int]] -> String
pShowShapes [] = ""
pShowShapes aL@(a:as) = "{" ++ pShowShape a ++ pShowShapesL as
    where pShowShapesL [] = "}"
          pShowShapesL (x:xs) = ", " ++ pShowShape x ++ pShowShapesL xs

pShowShape :: [Int] -> String
pShowShape [] = ""
pShowShape aL@(a:as) = "{" ++ show a ++ pShowShapeL as
    where pShowShapeL [] = "}"
          pShowShapeL (x:xs) = ", " ++ show x ++ pShowShapeL xs

pShowListIdentifiers :: [PName] -> String
pShowListIdentifiers [] = ""
pShowListIdentifiers (n:ns) = n ++ pShowListIdentifiersL ns
    where pShowListIdentifiersL [] = ""
          pShowListIdentifiersL nL@(n:ns) = ", " ++ pShowListIdentifiers nL

pShowDynamicDecl :: (Show a) => [([PName], PName, a)] -> (a -> String) -> String
pShowDynamicDecl [] _ = ""
pShowDynamicDecl (p:ps) showA = pShowDynamicDeclItem p showA ++ pShowDynamicDeclL ps showA
    where pShowDynamicDeclL [] _ = ""
          pShowDynamicDeclL xL@(x:xs) showA = ", " ++ pShowDynamicDecl xL showA

pShowDynamicDeclItem :: (Show a) => ([PName], PName, a) -> (a -> String) -> String
pShowDynamicDeclItem ([], v, attr) showA = 
    let sAttr = showA attr
    in  if sAttr == "" then v else v ++ " ( " ++ sAttr ++ " ) "
pShowDynamicDeclItem (qL@(q:qs), v, attr) showA = 
    let sAttr = showA attr
    in  if sAttr == "" then intercalate " " qL ++ v 
                       else intercalate " " qL ++ v ++ " ( " ++ sAttr ++ " ) "

pShowArrayDim :: [DimExpr] -> String
pShowArrayDim [] = ""
pShowArrayDim (c:cs) = show c ++ pShowArrayDimL cs
    where pShowArrayDimL [] = ""
          pShowArrayDimL (d:ds) = ", " ++ show d ++ pShowArrayDimL ds


