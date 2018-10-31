module G3 where

import GA

e1, e2, e3 :: Multivector Double Mother
e1 = toMulti [(Blade [Pos "e1" 1], (1))]
e2 = toMulti [(Blade [Pos "e2" 2], (1))]
e3 = toMulti [(Blade [Pos "e3" 3], (1))]
