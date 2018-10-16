module CGA where

{-
contains the definition for conformal geometric algebra
-}

import GA

data CGA = E1 | E2 | E3 | Ep | Em deriving (Eq, Ord, Show)

instance InnerProduct CGA where
    dot E1 E1 = 1
    dot E2 E2 = 1
    dot E3 E3 = 1
    dot Ep Ep = 1
    dot Em Em = -1
    dot _ _ = 0

zero, one :: MVector CGA
zero = toMulti []
one = toMulti [(Blade [], 1)]

scalar s = toMulti [(Blade [], s)]
vector e = toMulti [(Blade [e], 1)]

e1 :: MVector CGA
e1 = toMulti [(Blade [E1], 1)]
e2 = toMulti [(Blade [E2], 1)]
e3 = toMulti [(Blade [E3], 1)]
ep = toMulti [(Blade [Ep], 1)]
em = toMulti [(Blade [Em], 1)]

n0 = scale (1/2) (add ep em)
nInf = sub em ep

--project a R3 vector to CGA
proj a = (n0 `add` a `add` (scale (1/2) (a `mult` a `mult` nInf)))

--retract a = a `sub` n0 `sub` (scalar (1/2) `mult`)

--translate
trans a b = ((scalar 1) `sub` ((scalar (1/2)) `mult` a ` mult` nInf)) `mult` b `mult` ((scalar 1) `add` ((scalar (1/2)) `mult` a ` mult` nInf))

--rotation
--rot theta plane m =
