module CGA where

{-
contains the definition for conformal geometric algebra
-}

import GA

e1, e2, e3, em, ep, n0, nInf :: Multivector Sym Mother
e1 = toMulti [(Blade [Pos "e1" 1], (Lit 1))]
e2 = toMulti [(Blade [Pos "e2" 2], (Lit 1))]
e3 = toMulti [(Blade [Pos "e3" 3], (Lit 1))]
ep = toMulti [(Blade [Pos "e+" 4], (Lit 1))]
em = toMulti [(Blade [Neg "e-" 5], (Lit 1))]

scalar s = toMulti [(Blade [], s)]
vector e = toMulti [(Blade [e], 1)]

n0 = scale (Lit $ 1/2) (mAdd ep em)
nInf = mSub em ep

--project a R3 vector to CGA
proj a = (n0 `mAdd` a `mAdd` (scale (Lit $ 1/2) (a `mMul` a `mMul` nInf)))

--retract a = a `sub` n0 `sub` (scalar (1/2) `mult`)

--translate
--trans a b = ((scalar 1) `sub` ((scalar (1/2)) `mult` a ` mult` nInf)) `mult` b `mult` ((scalar 1) `add` ((scalar (1/2)) `mult` a ` mult` nInf))

--rotation
--rot theta plane m =
