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

getFromStmts :: (PArray -> Expr -> [Iter]) -> Map.Map PName PArray -> [Stmt] -> [Iter]
getFromStmts l_action _ [] = []
getFromStmts l_action l_arrayMap l_stmts@(a:as) = 
    let i1 = getFromStmt a 
        i2 = getFromStmts l_action l_arrayMap as 
    in  union i1 i2
    where getFromStmt (BRACES stmts) = getFromStmts l_action l_arrayMap stmts 
          getFromStmt (EXPR e) = getFromExpr e
          getFromStmt (DEXPR qs t es) = concat $ map getFromExpr es
          getFromStmt (IF e s1 s2) = 
              let iter1 = getFromExpr e 
                  iter2 = getFromStmt s1 
                  iter3 = getFromStmt s2 
              in  union iter1 (union iter2 iter3)
          getFromStmt (SWITCH e stmts) = 
              let iter1 = getFromExpr e 
                  iter2 = getFromStmts l_action l_arrayMap stmts
              in  union iter1 iter2
          getFromStmt (CASE v stmts) = getFromStmts l_action l_arrayMap stmts
          getFromStmt (DEFAULT stmts) = getFromStmts l_action l_arrayMap stmts
          getFromStmt (NOP) = []
          getFromStmt (BREAK) = []
          getFromStmt (DO e stmts) = 
              let iter1 = getFromExpr e 
                  iter2 = getFromStmts l_action l_arrayMap stmts
              in  union iter1 iter2
          getFromStmt (WHILE e stmts) = 
              let iter1 = getFromExpr e 
                  iter2 = getFromStmts l_action l_arrayMap stmts
              in  union iter1 iter2
          getFromStmt (FOR sL s) = 
              let iter1 = getFromStmt s 
                  iter2 = concat $ map (getFromStmts l_action l_arrayMap) sL
              in  union iter1 iter2
          getFromStmt (CONT) = []
          getFromStmt (RET e) = getFromExpr e
          getFromStmt (RETURN) = []
          getFromStmt (UNKNOWN s) = []
          getFromExpr (VAR q v) = []
          getFromExpr (PVAR q v dL) = 
              case Map.lookup v l_arrayMap of
                   Nothing -> []
                   Just arrayInUse -> l_action arrayInUse (PVAR q v dL)
          getFromExpr (BVAR v dim) = []
          getFromExpr (BExprVAR v e) = getFromExpr e
          getFromExpr (SVAR t e c f) = getFromExpr e
          getFromExpr (PSVAR t e c f) = getFromExpr e
          getFromExpr (Uno uop e) = getFromExpr e
          getFromExpr (PostUno uop e) = getFromExpr e
          getFromExpr (Duo bop e1 e2) = 
              let iter1 = getFromExpr e1 
                  iter2 = getFromExpr e2
              in  (union iter1 iter2)
          getFromExpr (PARENS e) = getFromExpr e
          getFromExpr _ = []

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
pIterLookup (v, dL) ((iterName, arrayInUse, dim):is) 
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
                l_arrayMacroName = l_arrayName ++ l_macro
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

pShowKernel :: String -> PKernel -> String
pShowKernel l_name l_kernel = "Pochoir_Kernel_" ++ show dim ++ "D(" ++ l_name ++ ", " ++
    pShowKernelParams (kParams l_kernel) ++ ")" ++ show (kStmt l_kernel) ++
    breakline ++ "Pochoir_Kernel_end" ++ breakline
        where dim = length (kParams l_kernel) - 1

-- AutoKernel is a de-sugared kernel
pShowAutoKernel :: String -> PKernel -> String
pShowAutoKernel l_name l_kernel = 
    let l_params = zipWith (++) (repeat "int ") (kParams l_kernel)
    in  "/* known! */ auto " ++ l_name ++ " = [&] (" ++ 
        pShowKernelParams l_params ++ ") {" ++ breakline ++
        show (kStmt l_kernel) ++
        breakline ++ "};" ++ breakline

pShowKernelParams :: [String] -> String
pShowKernelParams l_kernel_params = intercalate ", " l_kernel_params

pShowArrayGaps :: Int -> [PArray] -> String
pShowArrayGaps _ [] = ""
pShowArrayGaps l_rank l_array = breakline ++ "int " ++ 
        intercalate ", " (map (getArrayGaps (l_rank-1)) l_array) ++ ";"

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
        l_array = unionArrayIter l_iter
    --  in  pShowShadowArrayInUse l_array ++ pShowAutoKernel l_name l_kernel 
    in  pShowShadowArrayInUse l_sArrayInUse ++ pShowAutoKernel l_name l_kernel 

pShowMacroKernel :: PName -> [PArray] -> String -> PKernel -> String
pShowMacroKernel l_macro l_sArrayInUse l_name l_kernel =
    let l_iter = kIter l_kernel
        l_array = unionArrayIter l_iter
        -- shadowArrayInUse = pDefMacroArrayInUse l_array (kParams l_kernel)
        shadowArrayInUse = pDefMacroArrayInUse l_macro l_sArrayInUse (kParams l_kernel)
        -- unshadowArrayInUse = pUndefMacroArrayInUse l_array (kParams l_kernel)
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
        pShowIters l_iter ++ pShowArrayGaps l_rank l_array ++
        breakline ++ pShowStrides l_rank l_array ++ breakline ++
        "for (int " ++ l_t ++ " = t0; " ++ l_t ++ " < t1; ++" ++ l_t ++ ") { " ++ 
        pShowIterSet l_iter (kParams l_kernel)++
        breakline ++ pShowObaseForHeader l_rank l_iter (tail $ kParams l_kernel) ++
        breakline ++ pShowObaseStmt l_kernel ++ breakline ++ pShowObaseForTail l_rank ++
        pShowObaseTail l_rank ++ breakline ++ "};\n"

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
        breakline ++ pShowStrides l_rank l_array ++ breakline ++
        "for (int " ++ l_t ++ " = t0; " ++ l_t ++ " < t1; ++" ++ l_t ++ ") { " ++ 
        pShowPointerSet l_iter (kParams l_kernel)++
        breakline ++ pShowPointerForHeader l_rank l_iter (tail $ kParams l_kernel) ++
        breakline ++ pShowPointerStmt l_kernel ++ breakline ++ pShowObaseForTail l_rank ++
        pShowObaseTail l_rank ++ breakline ++ "};\n"

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
        breakline ++ pShowStrides l_rank l_array ++ breakline ++
        "for (int " ++ l_t ++ " = t0; " ++ l_t ++ " < t1; ++" ++ l_t ++ ") { " ++ 
        pShowOptPointerSet l_iter (kParams l_kernel)++
        breakline ++ pShowPointerForHeader l_rank l_iter (tail $ kParams l_kernel) ++
        breakline ++ pShowOptPointerStmt l_kernel ++ breakline ++ pShowObaseForTail l_rank ++
        pShowObaseTail l_rank ++ breakline ++ "};\n"

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
        breakline ++ pShowStrides l_rank l_array ++ breakline ++
        pShowRefMacro (kParams l_kernel) l_array ++
        "for (int " ++ l_t ++ " = t0; " ++ l_t ++ " < t1; ++" ++ l_t ++ ") { " ++ 
        breakline ++ pShowRawForHeader (tail $ kParams l_kernel) ++
        breakline ++ pShowCPointerStmt l_kernel ++ breakline ++ pShowObaseForTail l_rank ++
        pShowObaseTail l_rank ++ breakline ++ pShowRefUnMacro l_array ++ 
        "};\n"

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
        ") " ++ l_name ++ "_base[" ++ pGetTimeOffset l_toggle (DimVAR l_t) ++
        " * l_" ++ l_name ++ "_total_size + " ++ 
        (intercalate " + " $ zipWith pMul l_dims $ pStrideList l_name l_rank) ++ "]" ++
        breakline ++ breakline ++ pShowRefMacro l_kernelParams as

pStrideList :: PName -> Int -> [String]
pStrideList a 1 = ["l_stride_" ++ a ++ "_" ++ show 0]
pStrideList a r = ["l_stride_" ++ a ++ "_" ++ show (r-1)] ++ (pStrideList a $ r-1)

pMul :: String -> String -> String
pMul a b = "(" ++ a ++ ") * " ++ b

pShowArrayInfo :: [PArray] -> String
pShowArrayInfo [] = ""
pShowArrayInfo arrayInUse = foldr pShowArrayInfoItem "" arrayInUse
    where pShowArrayInfoItem l_arrayItem str =
            let l_type = aType l_arrayItem
                l_name = aName l_arrayItem
            in  str ++ breakline ++ show l_type ++ " * " ++ l_name ++ "_base"  ++ 
                " = " ++ l_name ++ ".data();" ++ breakline ++
                "const int " ++ "l_" ++ l_name ++ "_total_size = " ++ l_name ++
                ".total_size();" ++ breakline

pShowStrides :: Int -> [PArray] -> String
pShowStrides n [] = ""
pShowStrides n aL@(a:as) = "const int " ++ getStrides n aL ++ ";\n"
    where getStrides n aL@(a:as) = intercalate ", " $ concat $ map (getStride n) aL
          getStride 1 a = let r = 0 
                          in  ["l_stride_" ++ (aName a) ++ "_" ++ show r ++
                              " = " ++ (aName a) ++ ".stride(" ++ show r ++ ")"]
          getStride n a = let r = n-1
                          in  ["l_stride_" ++ (aName a) ++ "_" ++ show r ++
                              " = " ++ (aName a) ++ ".stride(" ++ show r ++ ")"] ++
                              getStride (n-1) a

pShowPointers :: [Iter] -> String
pShowPointers [] = ""
pShowPointers iL@(i:is) = foldr pShowPointer "" iL
    where pShowPointer (nameIter, arrayInUse, dL) str =
                str ++ breakline ++ (show $ aType arrayInUse) ++ " * " ++ nameIter ++ ";"

pShowPointerStmt :: PKernel -> String
pShowPointerStmt l_kernel = 
    let oldStmts = kStmt l_kernel
        l_iter = kIter l_kernel
        obaseStmts = transStmts oldStmts $ transPointer l_iter
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

transPointer :: [Iter] -> Expr -> Expr
transPointer l_iters (PVAR q v dL) =
    case pPointerLookup (v, dL) l_iters of
        Nothing -> PVAR q v dL
        Just (iterName, arrayInUse, des) -> 
            BVAR iterName de
                where de = simplifyDimExpr naive_de
                      naive_de = foldr plusCombDimExpr x $ zipWith mulDimExpr strideL $ tail $ excludeDimExpr dL des 
                      strideL = pGetArrayStrideList (aRank arrayInUse) (aName arrayInUse)
                      x = (DimINT 0)
transPointer l_iters e = e

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
pPointerLookup (v, dL) ((iterName, arrayInUse, dL'):is)
    | v == aName arrayInUse && head dL == head dL' = Just (iterName, arrayInUse, dL')
    | otherwise = pPointerLookup (v, dL) is

transIter :: [Iter] -> Expr -> Expr
transIter l_iters (PVAR q v dL) =
    case pIterLookup (v, dL) l_iters of
        Nothing -> PVAR q v dL
        Just iterName -> VAR q iterName
transIter l_iters e = e

pShowIterSet :: [Iter] -> [PName] -> String
pShowIterSet iL@(i:is) l_kernelParams = concat $ map pShowIterSetTerm iL
    where pShowIterSetTerm (name, array, dim) = 
            breakline ++ name ++ ".set(" ++ show (pShowTransDim dim l_kernelParams) ++ ");" 
            --
-- PName : list of kernel parameters
pShowOptPointerSet :: [Iter] -> [PName] -> String
pShowOptPointerSet [] _ = ""
pShowOptPointerSet iL@(i:is) l_kernelParams = 
    let baseIters = transIterN 0 $ getBaseIter l_kernelParams iL
    in  pShowPointers baseIters ++ (concat $ map pShowOptPointerSetTerm baseIters) ++ pShowNonBaseIters baseIters iL
        where pShowOptPointerSetTerm (iterName, array, dim) = 
                let l_arrayName = aName array
                    l_arrayBaseName = l_arrayName ++ "_base"
                    l_arrayTotalSize = "l_" ++ l_arrayName ++ "_total_size"
                    l_arrayStrideList = 
                        pGetArrayStrideList (length l_kernelParams - 1) l_arrayName
                    l_transDimList = tail $ pShowTransDim dim l_kernelParams
                    l_arraySpaceOffset = 
                        intercalate " + " $ zipWith pCombineDim l_transDimList l_arrayStrideList
                    l_arrayTimeOffset = (pGetTimeOffset (aToggle array) (head dim)) ++ 
                                        " * " ++ l_arrayTotalSize
                in  breakline ++ iterName ++ " = " ++ l_arrayBaseName ++ " + " ++ 
                    l_arrayTimeOffset ++ " + " ++ l_arraySpaceOffset ++ ";" 

pShowNonBaseIters :: [Iter] -> [Iter] -> String
pShowNonBaseIters _ [] = ""
pShowNonBaseIters bL iL@(i:is) = pShowShiftFromBase i bL ++ pShowNonBaseIters bL is
    where pShowShiftFromBase _ [] = "/* NO baseIter found! */"
          pShowShiftFromBase (iName, iArray, iDims) ((bName, bArray, bDims):bs)
            | iArray == bArray && head iDims == head bDims =
                let l_arrayName = aName iArray
                    l_arrayStrideList =
                        pGetArrayStrideList (length iDims - 1) l_arrayName
                    l_transDimList = map simplifyDimExpr $ zipWith excludeBaseDim (tail iDims) (tail bDims)
                    l_arraySpaceOffset = 
                        intercalate " + " $ zipWith pCombineDim l_transDimList l_arrayStrideList
                in  breakline ++ iName ++ " = "  ++ bName ++ " + " ++ l_arraySpaceOffset ++ ";"
            | otherwise = pShowShiftFromBase (iName, iArray, iDims) bs

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
pShowPointerSet iL@(i:is) l_kernelParams = concat $ map pShowPointerSetTerm iL
    where pShowPointerSetTerm (iterName, array, dim) = 
            let l_arrayName = aName array
                l_arrayBaseName = l_arrayName ++ "_base"
                l_arrayTotalSize = "l_" ++ l_arrayName ++ "_total_size"
                l_arrayStrideList = 
                    pGetArrayStrideList (length l_kernelParams - 1) l_arrayName
                l_transDimList = tail $ pShowTransDim dim l_kernelParams
                l_arraySpaceOffset = 
                    intercalate " + " $ zipWith pCombineDim l_transDimList l_arrayStrideList
                l_arrayTimeOffset = (pGetTimeOffset (aToggle array) (head dim)) ++ 
                                    " * " ++ l_arrayTotalSize
            in  breakline ++ iterName ++ " = " ++ l_arrayBaseName ++ " + " ++ 
                l_arrayTimeOffset ++ " + " ++ l_arraySpaceOffset ++ ";" 

pGetTimeOffset :: Int -> DimExpr -> String
pGetTimeOffset toggle tDim 
    | toggle == 2 = "((" ++ show tDim ++ ")" ++ " & 0x1" ++ ")"
    | toggle == 4 = "((" ++ show tDim ++ ")" ++ " & 0x11" ++ ")"
    | otherwise = "((" ++ show tDim ++ ") % " ++ show toggle ++ ")"

pCombineDim :: DimExpr -> String -> String
-- l_stride_pa_0 may NOT necessary be "1", 
-- plus that we have already set all strides to be of type "const int"
pCombineDim de stride = "(" ++ show de ++ ") * " ++ stride

pGetArrayStrideList :: Int -> PName -> [String]
pGetArrayStrideList 1 arr = ["l_stride_" ++ arr ++ "_" ++ show 0]
pGetArrayStrideList n arr = ("l_stride_" ++ arr ++ "_" ++ show (n-1)):(pGetArrayStrideList (n-1) arr)

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
        
pShowIters :: [Iter] -> String
pShowIters [] = ""
pShowIters ((l_name, l_array, l_dim):is) = 
    let l_type = aType l_array
        l_rank = aRank l_array
        l_toggle = aToggle l_array
        l_arrayName = aName l_array
    in breakline ++ "Pochoir_Iterator<" ++ show l_type ++ ", " ++ show l_rank ++ 
       ", " ++ show l_toggle ++ "> " ++
       l_name ++ "(" ++ l_arrayName ++ ");" ++ pShowIters is   

unionArrayIter :: [Iter] -> [PArray]
unionArrayIter [] = []
unionArrayIter iL@(i:is) = union (getArrayItem i) (unionArrayIter is)
    where getArrayItem (_, a, _) = [a]

getArrayIter :: [Iter] -> [PName]
getArrayIter [] = []
getArrayIter iL@(i:is) = (getArrayItem i) ++ (getArrayIter is)
    where getArrayItem (_, a, _) = [aName a]

pShowObaseForTail :: Int -> String
pShowObaseForTail n 
    | n == 0 = "/* end for (sub-trapezoid) */ "
    | otherwise = "} " ++ pShowObaseForTail (n-1)

pShowObaseTail :: Int -> String
pShowObaseTail n = 
    breakline ++ "/* Adjust sub-trapezoid! */" ++
    breakline ++ "for (int i = 0; i < " ++ show n ++ "; ++i) {" ++ 
    breakline ++ "\tl_grid.x0[i] += l_grid.dx0[i]; l_grid.x1[i] += l_grid.dx1[i];" ++
    breakline ++ "}" ++
    breakline ++ "} /* end for t */"

-- pL is the parameter list of original user supplied computing kernel
pShowObaseForHeader :: Int -> [Iter] -> [PName] -> String
pShowObaseForHeader _ _ [] = ""
pShowObaseForHeader 1 iL pL = 
                           breakline ++ pShowForHeader 0 (unionArrayIter iL) pL ++ 
                           pShowIterComma iL ++
                           breakline ++ intercalate (", " ++ breakline) 
                                        (map ((++) "++" . getIterName) iL) ++ ") {"
pShowObaseForHeader n iL pL = 
                           breakline ++ pShowForHeader (n-1) (unionArrayIter iL) pL ++ 
                           pShowIterComma iL ++
                           breakline ++ intercalate (", " ++ breakline) 
                                     (zipWith wrapIterInc
                                        (map (getArrayGap (n-1)) (getArrayIter iL))
                                        (map getIterName iL)) ++ 
                           ") {" ++ pShowObaseForHeader (n-1) iL pL
    where wrapIterInc gap iter = iter ++ ".inc(" ++ gap ++ ")"

-- pL is the parameter list of original user supplied computing kernel
pShowPointerForHeader :: Int -> [Iter] -> [PName] -> String
pShowPointerForHeader _ _ [] = ""
pShowPointerForHeader 1 iL pL = 
                           breakline ++ pShowPragma ++
                           breakline ++ pShowForHeader 0 (unionArrayIter iL) pL ++  
                           pShowIterComma iL ++
                           breakline ++ intercalate (", " ++ breakline) 
                                        (map ((++) "++" . getIterName) iL) ++ ") {"
--                                        (map ((flip (++) "+=1") . getIterName) iL) ++ ") {"

pShowPointerForHeader n iL pL = 
                           breakline ++ pShowForHeader (n-1) (unionArrayIter iL) pL ++ 
                           pShowIterComma iL ++
                           breakline ++ intercalate (", " ++ breakline) 
                                     (zipWith wrapIterInc
                                        (map (getArrayGap (n-1)) (getArrayIter iL))
                                        (map getIterName iL)) ++ 
                           ") {" ++ pShowPointerForHeader (n-1) iL pL
    where wrapIterInc gap iter = iter ++ " += " ++ gap 

pShowIterComma :: [Iter] -> String
pShowIterComma [] = ""
pShowIterComma iL@(i:is) = ", "

pShowForHeader :: Int -> [PArray] -> [PName] -> String
pShowForHeader i _ [] = ""
pShowForHeader i aL pL = 
    let len_pL = length pL
        idx = pL !! (len_pL - 1 - i)
        l_rank = show i
    in  adjustGap i aL ++ "for (int " ++ idx ++ 
        " = l_grid.x0[" ++ l_rank ++
        "]; " ++ idx ++ " < l_grid.x1[" ++ l_rank ++ "]; ++" ++ idx 
--        "]; " ++ idx ++ " < l_grid.x1[" ++ l_rank ++ "]; " ++ idx ++ "+=1"
                    
adjustGap :: Int -> [PArray] -> String
adjustGap i [] = ""
adjustGap i aL@(a:as) = 
    if i > 0 then pShowAdjustGap i aL
             else ""
    where pShowAdjustGap i [] = ""
          pShowAdjustGap i aL@(a:as) = concat $ map (pShowAdjustGapTerm i) aL
          pShowAdjustGapTerm i a = "gap_" ++ (aName a) ++ "_" ++ show i ++ " = " ++
                                   "l_stride_" ++ (aName a) ++ "_" ++ show i ++ 
                                   " + (l_grid.x0[" ++ show (i-1) ++ "] - l_grid.x1[" ++
                                   show (i-1) ++ "]) * l_stride_" ++ (aName a) ++ "_" ++ 
                                   show (i-1) ++ ";" ++ breakline
          
getIterName :: Iter -> String
getIterName (name, _, _) = name

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


