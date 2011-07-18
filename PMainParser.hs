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
import PParser2
import PUtils
import PData
import PShow
-- import Text.Show
import Data.List
import qualified Data.Map as Map

pParser :: GenParser Char ParserState String
pParser = do tokens0 <- many $ pToken
             eof
             return $ concat tokens0
             -- start a second pass!
--             setInput $ concat tokens0
--             tokens1 <- many pToken1
--             return $ concat tokens1

pToken :: GenParser Char ParserState String
pToken = 
        try pParseCPPComment
    <|> try pParseMacro
    <|> try pParsePochoirArray
    <|> try pParsePochoirStencil
    <|> try pParsePochoirShapeInfo
    <|> try pParsePochoirDomain
    <|> try pParsePochoirAutoKernel
    <|> try pParsePochoirAutoGuard
    <|> try pParsePochoirArrayMember
    <|> try pParsePochoirStencilMember
    <|> do ch <- anyChar
           return [ch]
    <?> "line"

pParseMacro :: GenParser Char ParserState String
pParseMacro = 
    do reserved "#define"
       l_name <- identifier
       pMacroValue l_name

pParsePochoirArray :: GenParser Char ParserState String
pParsePochoirArray =
    do reserved "Pochoir_Array"
       (l_type, l_rank) <- angles $ try pDeclStatic
       l_arrayDecl <- commaSep1 pDeclDynamic
       l_delim <- pDelim 
       updateState $ updatePArray $ transPArray (l_type, l_rank) l_arrayDecl
       return (breakline ++ "/* Known*/ Pochoir_Array <" ++ show l_type ++ 
               ", " ++ show l_rank ++ "> " ++ 
               pShowDynamicDecl l_arrayDecl pShowArrayDim ++ l_delim)

pParsePochoirStencil :: GenParser Char ParserState String
pParsePochoirStencil = 
    do reserved "Pochoir"
       l_rank <- angles exprDeclDim
       l_rawStencils <- commaSep1 pDeclPochoir
       l_delim <- pDelim
       l_state <- getState
       let l_stencils = map pSecond l_rawStencils
       updateState $ updatePStencil $ transPStencil l_rank l_stencils
       return (breakline ++ "/* Known */ Pochoir <" ++ show l_rank ++ 
               "> " ++ pShowDynamicDecl l_rawStencils (intercalate ", ") ++ 
               l_delim ++ breakline)

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
       updateState $ updatePShape (l_name, l_rank, l_len, l_toggle, l_slopes, l_timeShift, l_shapes)
       return (breakline ++ "/* Known */ Pochoir_Shape <" ++ show l_rank ++ "> " ++ l_name ++ " [" ++ show l_len ++ "] = " ++ pShowShapes l_shapes ++ ";\n" ++ breakline ++ "/* toggle: " ++ show l_toggle ++ "; slopes: " ++ show l_slopes ++ " */" ++ breakline)

pParsePochoirDomain :: GenParser Char ParserState String
pParsePochoirDomain =
    do reserved "Pochoir_Domain"
       l_rangeDecl <- commaSep1 pDeclDynamic
       semi
       updateState $ updatePRange $ transURange l_rangeDecl
       return (breakline ++ "Pochoir_Domain " ++ 
               pShowDynamicDecl l_rangeDecl pShowArrayDim ++ ";\n")

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
           -- return ("/* comment */")
           return ("/*" ++ str ++ "*/")
    <|> do try (string "//")
           str <- manyTill anyChar (try $ eol)
           -- return ("// comment\n")
           return ("//" ++ str ++ "\n")

pParsePochoirAutoKernel :: GenParser Char ParserState String
pParsePochoirAutoKernel =
    do reserved "auto"
       l_kernel_name <- identifier
       reservedOp "="
       symbol "[&]"
       l_kernel_params <- parens $ commaSep1 (reserved "int" >> identifier)
       reserved "{"
       exprStmts <- manyTill pStatement (try $ reserved "};")
       let l_kernel = PKernel { kName = l_kernel_name, kTimeShift = 0, 
                                kParams = l_kernel_params,
                                kStmt = exprStmts, kIter = [] }
       updateState $ updatePKernel l_kernel
       return (pShowAutoKernel l_kernel_name l_kernel) 

pParsePochoirAutoGuard :: GenParser Char ParserState String
pParsePochoirAutoGuard =
    do reserved "auto"
       l_guard_name <- identifier
       reservedOp "="
       symbol "[&]"
       l_guard_params <- parens $ commaSep1 (reserved "int" >> identifier)
       reservedOp "->" 
       reserved "bool" 
       reserved "{"
       exprStmts <- manyTill pStatement (try $ reserved "};")
       let l_guard = PGuard { gName = l_guard_name, gParams = l_guard_params,
                                gStmt = exprStmts, gIter = [] }
       updateState $ updatePGuard l_guard
       return (pShowAutoGuard l_guard_name l_guard) 

pMacroValue :: String -> GenParser Char ParserState String
pMacroValue l_name = 
                -- l_value <- liftM fromInteger $ try (natural)
              do l_value <- try (natural) >>= return . fromInteger
           -- Because Macro is usually just 1 line, we omit the state update 
                 updateState $ updatePMacro (l_name, l_value)
                 return ("#define " ++ l_name ++ " " ++ show (l_value) ++ "\n")
          <|> do l_value <- manyTill anyChar $ try eol
                 return ("#define " ++ l_name ++ " " ++ l_value ++ "\n")
          <?> "Macro Definition"

transPArray :: (PType, Int) -> [([PName], PName, [DimExpr])] -> [(PName, PArray)]
transPArray (l_type, l_rank) [] = []
transPArray (l_type, l_rank) (p:ps) =
    let l_name = pSecond p
        l_dims = pThird p
    in  (l_name, PArray {aName = l_name, aType = l_type, aRank = l_rank, aDims = l_dims, aMaxShift = 0, aToggle = 0, aRegBound = False}) : transPArray (l_type, l_rank) ps

transPStencil :: Int -> [PName] -> [(PName, PStencil)]
transPStencil l_rank [] = []
-- sToggle by default is two (2)
transPStencil l_rank (p:ps) = 
    (p, PStencil 
        {sName = p, sRank = l_rank, sToggle = 0, sUnroll = 1, 
         sTimeShift = 0, sArrayInUse = [], sShape = emptyShape, 
         sRegBound = False}) : transPStencil l_rank ps

transURange :: [([PName], PName, [DimExpr])] -> [(PName, PRange)]
transURange [] = []
transURange (p:ps) = (l_name, PRange {rName = l_name, rFirst = l_first, rLast = l_last, rStride = DimINT 1}) : transURange ps
    where l_name = pSecond p
          l_first = head $ pThird p
          l_last = head . tail $ pThird p


