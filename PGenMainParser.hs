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

module PGenMainParser where

import Text.ParserCombinators.Parsec

import Control.Monad

import PGenBasicParser
import PGenMainParser2
import PGenUtils
import PGenData
import PGenShow
-- import Text.Show
import Data.List
import qualified Data.Map as Map

-- first pass will be gathering the info
-- second pass will do the real transformation/code generation
pParser :: GenParser Char ParserState String
pParser = do tokens0 <- many $ pToken
             eof
             -- return $ concat tokens0
             -- start a second pass!
             -- setInput $ concat tokens0
             -- tokens1 <- many pToken1
             l_state <- getState
             let l_mode = pMode l_state
             let l_colorNum = pColorNum l_state
             -- let l_stencils = Map.elems $ pGenPlan l_state
             let l_colorVectors = pColorVectors l_state
             let l_guardFuncs = Map.elems $ pGuardFunc l_state
             -- let l_reg_GTs = sRegTileKernel l_stencil
             case Map.lookup l_colorNum $ pGenPlan l_state of
                 Nothing -> 
                    return ("pGenPlan (with colorNum = " ++ (show l_colorNum) ++ 
                           ") Not Found!")
                 Just l_stencil -> 
                    do let l_arrayInUse = sArrayInUse l_stencil
                       let l_regBound = foldr (||) False $ map (getArrayRegBound l_state) l_arrayInUse
                       let l_stencil' = l_stencil { sRegBound = l_regBound }
                       let l_output = pCodeGen l_mode l_colorVectors l_guardFuncs l_stencil'
                       return $ fst l_output

pToken :: GenParser Char ParserState String
pToken = 
        try pParseCPPComment
    <|> try pParsePochoirArray
    <|> try pParsePochoirStencil
    <|> try pParsePochoirShapeInfo
    <|> try pParsePochoirDomain
    <|> try pParsePochoirTile
    <|> try pParsePochoirKernel
    <|> try pParsePochoirAutoKernelFunc
    <|> try pParsePochoirGuard
    <|> try pParsePochoirAutoGuardFunc
    <|> try pParsePochoirArrayMember
    <|> try pParsePochoirStencilMember
    <|> do ch <- anyChar
           -- return [ch]
           return [] 
    <?> "line"

pParsePochoirArray :: GenParser Char ParserState String
pParsePochoirArray =
    do reserved "Pochoir_Array"
       (l_type, l_rank) <- angles $ try pDeclStatic
       l_arrayDecl <- commaSep1 pDeclDynamic
       l_delim <- pDelim 
       updateState $ updatePArray $ transPArray (l_type, l_rank) l_arrayDecl
       return ""

pParsePochoirStencil :: GenParser Char ParserState String
pParsePochoirStencil = 
    do reserved "Pochoir"
       l_rank <- angles exprDeclDim
       l_rawStencils <- commaSep1 pDeclPochoir
       l_delim <- pDelim
       l_state <- getState
       let l_stencils = map pSecond l_rawStencils
       updateState $ updatePStencil $ transPStencil l_rank l_stencils
       return ""

pParsePochoirKernel :: GenParser Char ParserState String
pParsePochoirKernel =
    do reserved "Pochoir_Kernel"
       l_rank <- angles pDeclStaticNum
       l_name <- identifier
       (l_kernelFunc, l_shape) <- parens pPochoirKernelParams
       semi
       l_state <- getState
       case Map.lookup l_shape $ pShape l_state of
            Nothing -> return (breakline ++ "Pochoir_Kernel <" ++ 
                               show l_rank ++ "> " ++ l_name ++ "(" ++ 
                               l_shape ++ ", " ++ l_kernelFunc ++ 
                               "); /* UNKNOWN Shape */" ++ breakline)
            Just l_pShape ->
                case Map.lookup l_kernelFunc $ pKernelFunc l_state of
                     Nothing -> return ""
                     Just l_pKernelFunc -> 
                          do let l_kernel = PKernel { kName = l_name, kRank = l_rank, kShape = l_pShape, kFunc = l_pKernelFunc { kfShape = l_pShape, kfName = l_name }, kIndex = [], kComment = "" }
                             updateState $ updatePKernel l_kernel
                             return ""

pParsePochoirTile :: GenParser Char ParserState String
pParsePochoirTile =
    do reserved "Pochoir_Kernel"
       l_rank <- angles pDeclStaticNum
       l_name <- identifier
       l_sizes <- many1 $ brackets pDeclStaticNum
       reservedOp "="
       l_tile_kernel <- pParseTileKernel
       semi
       let l_tile = PTile { tName = l_name, tRank = l_rank, tSize = l_sizes, tKernel = l_tile_kernel, tComment = "", tOp = PNULL, tOrigGuard = emptyGuard, tOrder = 0, tColor = emptyColor }
       let l_kernels = getTileKernels l_tile
       updateState $ updatePTile l_tile
       return ""

pParsePochoirGuard :: GenParser Char ParserState String
pParsePochoirGuard =
    do reserved "Pochoir_Guard"
       l_rank <- angles pDeclStaticNum
       l_name <- identifier
       l_guardFunc <- parens identifier
       semi
       l_state <- getState
       case Map.lookup l_guardFunc $ pGuardFunc l_state of
            Nothing -> return ""
            Just l_pGuardFunc -> 
                 do let l_guard = emptyGuard { gName = l_name, gRank = l_rank, gFunc = l_pGuardFunc { gfName = l_name } }
                    updateState $ updatePGuard l_guard
                    return ""

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
                              shapeSlopes = l_slopes, 
                              shapeTimeShift = l_timeShift, 
                              shape = l_shapes, shapeComment = ""}
       updateState $ updatePShape l_pShape
       return ""

pParsePochoirDomain :: GenParser Char ParserState String
pParsePochoirDomain =
    do reserved "Pochoir_Domain"
       l_rangeDecl <- commaSep1 pDeclDynamic
       semi
       updateState $ updatePRange $ transURange l_rangeDecl
       return ""

pParsePochoirArrayMember :: GenParser Char ParserState String
pParsePochoirArrayMember =
    do l_id <- try (pIdentifier)
       l_state <- getState
       try $ ppArray l_id l_state

pParsePochoirStencilMember :: GenParser Char ParserState String
pParsePochoirStencilMember =
    do l_id <- try (pIdentifier)
       l_state <- getState
       try $ ppStencil l_id l_state

pParseCPPComment :: GenParser Char ParserState String
pParseCPPComment =
        do try (string "/*")
           str <- manyTill anyChar (try $ string "*/")
           return ""
    <|> do try (string "//")
           str <- manyTill anyChar (try $ eol)
           return ""

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
                                        kfTileOp = PNULL, kfGuardFunc = emptyGuardFunc,
                                        kfTileOrder = 0, kfComment = "" }
       updateState $ updatePKernelFunc l_kernelFunc
       return ""

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
       return ""

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
         sRegBound = False, sRegTileKernel = [], 
         sComment = ""}) : transPStencil l_rank ps

transURange :: [([PName], PName, [DimExpr])] -> [(PName, PRange)]
transURange [] = []
transURange (p:ps) = (l_name, PRange {rName = l_name, rFirst = l_first, rLast = l_last, rStride = DimINT 1, rComment = ""}) : transURange ps
    where l_name = pSecond p
          l_first = head $ pThird p
          l_last = head . tail $ pThird p

    
