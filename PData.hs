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

module PData where

import Text.ParserCombinators.Parsec
import Control.Monad

import Data.Char
import Data.List
import qualified Data.Map as Map

type PName = String
type PValue = Int
type Uop = String 
type Bop = String
breakline :: String
breakline = "\n\t"

stripWhite :: String -> String
stripWhite l = dropWhile isSpace l

-- We use newtype just because we want to manually derive a Show instance of PShift
newtype PShift = PShift Int deriving (Eq, Ord)
data RegionT = Periodic | Nonperiodic | UnknownRegionT deriving Show
data PBasicType = PChar | PShort | PLong | PSigned | PUnsigned | PVoid | PInt | PDouble | PFloat | PBool | PUserType deriving Eq
data PRWMode = PRead | PWrite | POther deriving (Eq, Show)
data PType = PType {
    basicType :: PBasicType,
    typeName :: String
} deriving Eq
data PState = PochoirBegin | PochoirEnd | PochoirMacro | PochoirDeclArray | PochoirDeclRange | PochoirError | Unrelated deriving (Show, Eq)
data PMode = PHelp | PDefault | PDebug | PCaching | PCPointer | POptPointer | PPointer | PMacroShadow | PNoPP deriving Eq
data PMacro = PMacro {
    mName :: PName,
    mValue :: PValue
} deriving Show
data PArray = PArray {
    aName :: PName,
    aType :: PType,
    aRank :: Int,
    aMaxShift :: Int,
    aToggle :: Int,
    aDims :: [DimExpr],
    aRegBound :: Bool
} deriving (Show, Eq)
data PStencil = PStencil {
    sName :: PName,
    sRank :: Int,
    sToggle :: Int,
    sUnroll :: Int,
    sArrayInUse :: [PArray],
    sShape :: PShape,
    sRegBound :: Bool
} deriving Show
data PShape = PShape {
    shapeName :: PName,
    shapeRank :: Int,
    shapeLen :: Int,
    shapeToggle :: Int,
    shapeSlopes :: [Int],
    shape :: [[Int]]
} deriving Show

emptyShape :: PShape
emptyShape = PShape {shapeName = "", shapeRank = 0, shapeLen = 0, shapeToggle = 0, shapeSlopes = [], shape = []}

data PRange = PRange {
    rName :: PName,
    rFirst :: DimExpr,
    rLast :: DimExpr,
    rStride :: DimExpr 
} deriving Show 

-- (NameOfIter, arrayInUse, correspondingDimExpr, RWMode)
type Iter = (String, PArray, [DimExpr], PRWMode)

data PKernel = PKernel {
    kName :: PName,
    kParams :: [PName],
    kStmt :: [Stmt],
    kIter :: [Iter]
} deriving Show

emptyKernel :: PKernel
emptyKernel = PKernel {kName = "", kParams = [], kStmt = [], kIter = []}

data ParserState = ParserState {
    pMode  :: PMode,
    pState :: PState, 
    pMacro :: Map.Map PName PValue, 
    pArray :: Map.Map PName PArray,
    pStencil :: Map.Map PName PStencil,
    pRange :: Map.Map PName PRange,
    pShape :: Map.Map PName PShape,
    pKernel :: Map.Map PName PKernel
} deriving Show

data Expr = VAR String String 
          -- PVAR : (*&)p(a, b, ...)
          | PVAR String String [DimExpr]
          -- BVAR : v[a]
          | BVAR String DimExpr
          -- BExprVAR : v[a+b+c]
          | BExprVAR String Expr
          -- SVAR : type(p(a, b)).field
          | SVAR PType Expr String String
          -- PSVAR : (type*)(p(a, b)).field
          | PSVAR PType Expr String String
          -- Uno is prefix unary operator
          | Uno Uop Expr
          -- PostUno is postfix unary operator
          | PostUno Uop Expr
          | Duo Bop Expr Expr
          | INT Int
          | FLOAT Double
          | BOOL String
          | PARENS Expr
          deriving Eq

data Stmt = BRACES [Stmt]
          | EXPR Expr
          | DEXPR [PName] PType [Expr]
          | IF Expr Stmt Stmt
          | SWITCH Expr [Stmt]
          | CASE PValue [Stmt]
          | DEFAULT [Stmt]
          | NOP
          | BREAK
          | DO Expr [Stmt]
          | WHILE Expr [Stmt]
          | FOR [[Stmt]] Stmt
          | CONT
          | RET Expr
          | RETURN
          | UNKNOWN String
          deriving Eq

data DimExpr = DimVAR String 
          | DimDuo Bop DimExpr DimExpr
          | DimParen DimExpr
          | DimINT Int
          deriving Eq

instance Show PMode where
    show PHelp = "-help" 
    show PDefault = "-default" 
    show PDebug = "-debug" 
    show PCaching = "-split-caching" 
    show PCPointer = "-split-c-pointer" 
    show POptPointer = "-split-opt-pointer" 
    show PPointer = "-split-pointer" 
    show PMacroShadow = "-split-macro-shadow" 
    show PNoPP = "-No-Preprocessing"

instance Show PType where
    show ptype = typeName ptype

instance Show Expr where
    show (VAR q str) = q ++ str
    show (PVAR q a xList@(x:xs)) = q ++ a ++ "(" ++ showList xList "" ++ ")"
    show (BVAR a x) = a ++ "[" ++ show x ++ "]"
    show (BExprVAR a e) = a ++ "[" ++ show e ++ "]"
    show (SVAR t e c f) = show t ++ "(" ++ show e ++ ")" ++ c ++ f
    show (PSVAR t e c f) = "(" ++ show t ++ ")(" ++ show e ++ ")" ++ c ++ f
    show (Uno uop expr) = uop ++ show expr 
    show (PostUno uop expr) = show expr ++ uop
    show (Duo bop lexpr rexpr) = show lexpr ++ " " ++ bop ++ " " ++ show rexpr
    show (PARENS expr) = "(" ++ show expr ++ ")"
    show (INT n) = show n
    show (FLOAT n) = show n
    show (BOOL b) = b
    showList [] = showString ""
    showList (x:xs) = showString breakline . shows x . showChar ';' . showList xs

instance Show DimExpr where
    show (DimVAR str) = str
    show (DimDuo bop lexpr rexpr) = show lexpr ++ " " ++ bop ++ " " ++ show rexpr
    show (DimParen e) = "(" ++ show e ++ ")"
    show (DimINT n) = show n
    showList (x:xs) = shows x . showl xs
                      where showl [] = showString "" 
                            showl (x:xs) = showString ", " . shows x . showl xs

instance Show Stmt where
    show (BRACES tL@(t:ts)) = "{" ++ showList tL "" ++ breakline ++ "}"
    show (EXPR expr) = show expr ++ ";"
    show (DEXPR qs declType es) = intercalate " " qs ++ " " ++ 
        (intercalate ", " $ map show es) ++ ";"
    show (IF expr l_stmt NOP) = "if " ++ show expr ++ 
        breakline ++ show l_stmt 
    show (IF expr l_stmt r_stmt) = "if (" ++ show expr ++ ") {" ++ 
        breakline ++ show l_stmt ++ 
        breakline ++ "} else {" ++ show r_stmt ++ "}"
    show (SWITCH expr tL@(t:ts)) = "switch " ++ show expr ++ "{" ++
        showList tL "" ++ breakline ++ "} /* end of switch */" ++ breakline
    show (CASE l_value tL@(t:ts)) = "case " ++ show l_value ++ " : " ++
        showList tL ""
    show (DEFAULT tL@(t:ts)) = "default :" ++ showList tL "" 
    show (NOP) = ""
    show (BREAK) = "break;"
    show (DO expr stmts) = "do {" ++ 
                            breakline ++ showList stmts "" ++ breakline ++ 
                            " }while " ++ show expr ++ ";" ++ breakline
    show (WHILE expr stmts) = "while (" ++ show expr ++ ") {" ++
        showList stmts "" ++ breakline ++ "} /* end of while */" ++ breakline
    show (FOR ttL@(t:ts) l_stmt) = "for " ++ showForListList ttL ++
        breakline ++ show l_stmt
            where showForListList (t:ts) = "(" ++ showForList t ++ showForListListL ts
                  showForListListL [] = ")"
                  showForListListL (t:ts) = "; " ++ showForList t ++ showForListListL ts
                  showForList (t:ts) = showForExpr t ++ showForListL ts
                  showForListL [] = ""
                  showForListL (t:ts) = ", " ++ showForExpr t ++ showForListL ts
                  showForExpr NOP = ""
                  showForExpr (EXPR expr) = show expr 
                  showForExpr (DEXPR qs declType expr) = intercalate " " qs ++ 
                                                         " " ++ show expr
    show (CONT) = "continue;" ++ breakline
    show (RET e) = "return (" ++ show e ++ ");" ++ breakline
    show (RETURN) = "return;" ++ breakline
    show (UNKNOWN stmt) = stmt ++ "\t/* Unrecognized! */" 
    showList [] = showString ""
    showList (x:xs) = shows x . showString breakline . showList xs
    -- showList (x:xs) = showString breakline . shows x . showList xs

instance Show PShift where
    show (PShift n) = show n
    showList [] = showString ""
    showList (x:xs) = showChar '{' . shows x . showl xs
                        where showl [] = showChar '}'
                              showl (x:xs) = showString ", " . shows x . showl xs


