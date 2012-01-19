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

module PMainParser where

import Text.ParserCombinators.Parsec

import Control.Monad

import PBasicParser
import PMainParser2
import PUtils
import PData
import PShow
-- import Text.Show
import Data.List
import qualified Data.Map as Map

-- first pass will be gathering the info
-- second pass will do the real transformation
pParser :: GenParser Char ParserState (String, String)
pParser = do tokens0 <- many $ pToken
             eof
             -- starting a second pass!
             -- the type of tokens0 is [(String, String)]
             setInput $ concatMap fst tokens0
             tokens1 <- many pToken1
             let out_info = concat tokens1
             let pochoir_info = concatMap snd tokens0
             return (out_info, pochoir_info)

pToken :: GenParser Char ParserState (String, String)
pToken = 
        do l_ret <- try pParseCPPComment
           return (l_ret, "")
    <|> do l_ret <- try pIncludeFile
           return (l_ret, l_ret)
    <|> do l_ret <- try pParsePochoirArray
           return (l_ret, l_ret)
    <|> do l_ret <- try pParsePochoirStencil
           return (l_ret, l_ret)
    <|> do l_ret <- try pParsePochoirShapeInfo
           return (l_ret, l_ret)
    <|> do l_ret <- try pParsePochoirDomain
           return (l_ret, l_ret)
    <|> do l_ret <- try pParsePochoirTile
           return (l_ret, l_ret)
    <|> do try pParsePochoirKernel
    <|> do l_ret <- try pParsePochoirAutoKernelFunc
           return (l_ret, l_ret)
    <|> do try pParsePochoirGuard
    <|> do l_ret <- try pParsePochoirAutoGuardFunc
           return (l_ret, l_ret)
    <|> do try pParsePochoirArrayMember
    <|> do try pParsePochoirStencilMember
    <|> do ch <- anyChar
           return ([ch], "")
    <?> "line"

pIncludeFile :: GenParser Char ParserState String
pIncludeFile = 
    do reserved "#include"
       try (angles identifier >>= \l_name -> return (breakline ++ "#include <" ++ l_name ++ ">" ++ breakline)) 
            <|> try (angles (identifier >>= \x -> symbol ".h" >> return (x ++ ".h")) >>= \l_name -> return (breakline ++ "#include <" ++ l_name ++ ">" ++ breakline)) 
            <|> try (symbol "\"" >> identifier >>= \x -> symbol ".h" >> return (x ++ ".h") >>= \x -> symbol "\"" >> return x >>= \l_name -> return (breakline ++ "#include \"" ++ l_name ++ "\"" ++ breakline)) 
            <|> try (symbol "\"" >> identifier >>= \x -> symbol "\"" >> return x >>= \l_name -> return (breakline ++ "#include \"" ++ l_name ++ "\"" ++ breakline))

pParsePochoirArray :: GenParser Char ParserState String
pParsePochoirArray =
    do reserved "Pochoir_Array"
       (l_type, l_rank) <- angles $ try pDeclStatic
       l_arrayDecl <- commaSep1 pDeclDynamic
       l_delim <- pDelim 
       updateState $ updatePArray $ transPArray (l_type, l_rank) l_arrayDecl
       let l_ret = breakline ++ "/* Known*/ Pochoir_Array <" ++ 
                   show l_type ++ ", " ++ show l_rank ++ "> " ++ 
                   pShowDynamicDecl l_arrayDecl pShowArrayDim ++ l_delim
       return l_ret

pParsePochoirStencil :: GenParser Char ParserState String
pParsePochoirStencil = 
    do reserved "Pochoir"
       l_rank <- angles exprDeclDim
       l_rawStencils <- commaSep1 pDeclPochoir
       l_delim <- pDelim
       l_state <- getState
       let l_stencils = map pSecond l_rawStencils
       updateState $ updatePStencil $ transPStencil l_rank l_stencils
       let l_ret = breakline ++ "/* Known */ Pochoir <" ++ show l_rank ++ 
               "> " ++ pShowDynamicDecl l_rawStencils (showString "") ++ 
               l_delim ++ breakline
       return l_ret

pParsePochoirKernel :: GenParser Char ParserState (String, String)
pParsePochoirKernel =
    do reserved "Pochoir_Kernel"
       l_rank <- angles pDeclStaticNum
       l_name <- identifier
       (l_kernelFunc, l_shape) <- parens pPochoirKernelParams
       semi
       l_state <- getState
       case Map.lookup l_shape $ pShape l_state of
            Nothing -> let l_ret = breakline ++ "Pochoir_Kernel <" ++ 
                               show l_rank ++ "> " ++ l_name ++ "(" ++ 
                               l_shape ++ ", " ++ l_kernelFunc ++ 
                               "); /* UNKNOWN Shape */" ++ breakline
                       in  return (l_ret, "")
            Just l_pShape ->
                case Map.lookup l_kernelFunc $ pKernelFunc l_state of
                     Nothing -> let l_ret = 
                                        breakline ++ "Pochoir_Kernel <" ++ 
                                        show l_rank ++ "> " ++ l_name ++ 
                                        "(" ++ l_shape ++ ", " ++ 
                                        l_kernelFunc ++ 
                                        "); /* UNKNOWN Kernel Func */ " ++ 
                                        breakline
                                in  return (l_ret, "")
                     Just l_pKernelFunc -> 
                          do let l_kernel = PKernel { kName = l_name, kRank = l_rank, kShape = l_pShape, kFunc = l_pKernelFunc { kfShape = l_pShape, kfName = l_name }, kIndex = [], kComment = "" }
                             updateState $ updatePKernel l_kernel
                             let l_ret = 
                                     breakline ++ "Pochoir_Kernel <" ++ 
                                     show l_rank ++ "> " ++ l_name ++ 
                                     "(" ++ l_kernelFunc ++ ", " ++ 
                                     l_shape ++ 
                                     "); /* KNOWN!!! */" ++ breakline
                             return (l_ret, l_ret)

pParsePochoirTile :: GenParser Char ParserState String
pParsePochoirTile =
    do reserved "Pochoir_Kernel"
       l_rank <- angles pDeclStaticNum
       l_name <- identifier
       l_sizes <- many1 $ brackets pDeclStaticNum
       reservedOp "="
       l_tile_kernel <- pParseTileKernel
       semi
       l_state <- getState
       -- We should defer the order to register_tile_x functions
       let l_tile = PTile { tName = l_name, tRank = l_rank, tSize = l_sizes, tKernel = l_tile_kernel, tComment = "", tOp = PNOP, tOrigGuard = emptyGuard, tOrder = 0 }
       let l_kernels = getTileKernels l_tile
       updateState $ updatePTile l_tile
       let l_ret = 
            breakline ++ "Pochoir_Kernel <" ++ show l_rank ++ "> " ++ 
            l_name ++ pShowArrayDims l_sizes ++ " = " ++ 
            show l_tile_kernel ++ ";" ++ " /* Known! */" ++ breakline ++ 
            "/* " ++ (intercalate "; " $ 
                        zipWith (++) (map kName l_kernels) 
                                     (map (show . kIndex) l_kernels)) ++ 
            " */" ++ breakline
       return l_ret

pParsePochoirGuard :: GenParser Char ParserState (String, String)
pParsePochoirGuard =
    do reserved "Pochoir_Guard"
       l_rank <- angles pDeclStaticNum
       l_name <- identifier
       l_guardFunc <- parens identifier
       semi
       l_state <- getState
       case Map.lookup l_guardFunc $ pGuardFunc l_state of
            Nothing -> let l_ret = breakline ++ "Pochoir_Guard <" ++ 
                               show l_rank ++ "> " ++ l_name ++ "(" ++ 
                               l_guardFunc ++ 
                               "); /* UNKNOWN Guard Func */" ++ breakline
                       in  return (l_ret, "")
            Just l_pGuardFunc -> 
                 do let l_guard = emptyGuard { gName = l_name, gRank = l_rank, gFunc = l_pGuardFunc { gfName = l_name } }
                    updateState $ updatePGuard l_guard
                    let l_ret = breakline ++ "Pochoir_Guard <" ++ 
                            show l_rank ++ "> " ++ l_name ++ 
                            "(" ++ l_guardFunc ++ "); /* Known!!! */" ++ 
                            breakline
                    return (l_ret, l_ret)

pParsePochoirShapeInfo :: GenParser Char ParserState String
pParsePochoirShapeInfo = 
    do reserved "Pochoir_Shape"
       l_rank <- angles pDeclStaticNum
       l_name <- identifier
       brackets $ option 0 pDeclStaticNum
       reservedOp "="
       l_shapes <- braces (commaSep1 ppShape)
       semi
       let l_len = length l_shapes
       let l_toggle = getToggleFromShape l_shapes
       let l_slopes = getSlopesFromShape (l_toggle-1) l_shapes 
       let l_timeShift = getTimeShiftFromShape l_shapes
       let l_pShape = PShape {shapeName = l_name, shapeRank = l_rank, 
                              shapeLen = l_len, shapeToggle = l_toggle, 
                              shapeSlopes = l_slopes, shapeTimeShift = l_timeShift, 
                              shape = l_shapes, shapeComment = ""}
       updateState $ updatePShape l_pShape
       return (show l_pShape)

pParsePochoirDomain :: GenParser Char ParserState String
pParsePochoirDomain =
    do reserved "Pochoir_Domain"
       l_rangeDecl <- commaSep1 pDeclDynamic
       semi
       updateState $ updatePRange $ transURange l_rangeDecl
       return (breakline ++ "Pochoir_Domain " ++ 
               pShowDynamicDecl l_rangeDecl pShowArrayDim ++ ";\n")

pParsePochoirArrayMember :: GenParser Char ParserState (String, String)
pParsePochoirArrayMember =
    do l_id <- try (pIdentifier)
       l_state <- getState
       try $ ppArray l_id l_state

pParsePochoirStencilMember :: GenParser Char ParserState (String, String)
pParsePochoirStencilMember =
    do l_id <- try (pIdentifier)
       l_state <- getState
       try $ ppStencil l_id l_state

pParseCPPComment :: GenParser Char ParserState String
pParseCPPComment =
        do try (string "/*")
           str <- manyTill anyChar (try $ string "*/")
           -- return ("/* comment */")
           return ("/*" ++ str ++ "*/")
    <|> do try (string "//")
           str <- manyTill anyChar (try $ eol)
           -- return ("// comment\n")
           return ("//" ++ str ++ "\n")

pParsePochoirAutoKernelFunc :: GenParser Char ParserState String
pParsePochoirAutoKernelFunc =
    do reserved "auto"
       l_kernel_name <- identifier
       reservedOp "="
       symbol "[&]"
       l_kernel_params <- parens $ commaSep1 (reserved "int" >> identifier)
       reserved "{"
       exprStmts <- manyTill pStatement (try $ reserved "};")
       let l_len_kernel_params = length l_kernel_params
       -- The transformation of input kernel parameters (t, i, j, ...) to 
       -- (t, i2, i1, ...) is not safe because besides the PVAR ..., 
       -- the input kernel parameters t, i, j, ... can also be used in other 
       -- calculation spread through the kernel functions.
       let l_rev_exprStmts = transStmts exprStmts (transKernelParams l_kernel_params 0 (l_len_kernel_params-1))
       let l_rev_kernel_params = getKernelParams l_len_kernel_params
       let l_kernelFunc = PKernelFunc { kfName = l_kernel_name, 
                                        kfParams = l_kernel_params,
                                        kfStmt = exprStmts, 
                                        kfStmtSize = length exprStmts,
                                        kfIter = [], kfShape = emptyShape, 
                                        kfTileOp = PNOP, kfGuardFunc = emptyGuardFunc,
                                        kfTileOrder = 0, kfComment = "" }
       updateState $ updatePKernelFunc l_kernelFunc
       let l_ret = pShowAutoKernelFunc l_kernel_name l_kernelFunc
       return l_ret

pParsePochoirAutoGuardFunc :: GenParser Char ParserState String
pParsePochoirAutoGuardFunc =
    do reserved "auto"
       l_guard_name <- identifier
       reservedOp "="
       symbol "[&]"
       l_guard_params <- parens $ commaSep1 (reserved "int" >> identifier)
       reservedOp "->" 
       reserved "bool" 
       reserved "{"
       exprStmts <- manyTill pStatement (try $ reserved "};")
       let l_guardFunc = PGuardFunc { gfName = l_guard_name, gfParams = l_guard_params,
                                      gfStmt = exprStmts, 
                                      gfStmtSize = length exprStmts, 
                                      gfIter = [], gfComment = "" }
       updateState $ updatePGuardFunc l_guardFunc
       let l_ret = pShowAutoGuardFunc l_guard_name l_guardFunc
       return l_ret

transPArray :: (PType, Int) -> [([PName], PName, [DimExpr])] -> [(PName, PArray)]
transPArray (l_type, l_rank) [] = []
transPArray (l_type, l_rank) (p:ps) =
    let l_name = pSecond p
        l_dims = pThird p
    in  (l_name, PArray {aName = l_name, aType = l_type, aRank = l_rank, aDims = l_dims, aMaxShift = 0, aToggle = 0, aRegBound = False, aComment = ""}) : transPArray (l_type, l_rank) ps

transPStencil :: Int -> [PName] -> [(PName, PStencil)]
transPStencil l_rank [] = []
-- sToggle by default is two (2)
transPStencil l_rank (p:ps) = 
    (p, PStencil 
        {sName = p, sRank = l_rank, sToggle = 0, 
         sTimeShift = 0, sArrayInUse = [], sShape = emptyShape, 
         sRegBound = False, sRegStaggerKernel = [], 
         sRegTileKernel = [], sRegInclusiveTileKernel = [], 
         sRegTinyInclusiveTileKernel = [], 
         sComment = ""}) : transPStencil l_rank ps

transURange :: [([PName], PName, [DimExpr])] -> [(PName, PRange)]
transURange [] = []
transURange (p:ps) = (l_name, PRange {rName = l_name, rFirst = l_first, rLast = l_last, rStride = DimINT 1, rComment = ""}) : transURange ps
    where l_name = pSecond p
          l_first = pHead $ pThird p
          l_last = pHead . tail $ pThird p

    
