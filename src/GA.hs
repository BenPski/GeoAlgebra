{-# LANGUAGE MultiParamTypeClasses #-}
module GA where

import qualified Data.Map.Strict as M
import Data.List (intercalate)

{-

Want a nicer representation for multivectors
this way canonicalizing the form doesn't require doing psedo raising of the blades to multivectors

Notably everything is much simpler when restriciting the generating vectors to being orthonomal (a subalgebra of R(inf,inf), the mother algebra)
This means that arbitrary bases do not work (e.g., directly defined by null vectors), but they can be operated on like they are by first defining the orthonomral basis and then defining and using the null vectors

-}

--defined on a scalar and the generating vectors
newtype Multivector s a = Multivector (M.Map (Blade a) s) deriving (Eq, Ord)

--split into a list of terms
toTerms :: Multivector s a -> [Multivector s a]
toTerms (Multivector a) = fmap (\(b,v) -> Multivector $ M.singleton b v) (M.assocs a)

fromTerms :: (Eq s, Ring s, InnerProduct a s, Ord a) => [Multivector s a] -> Multivector s a
fromTerms ms = makeMulti $ concatMap (\(Multivector a) -> M.assocs a) ms

instance (Show a, Show s, Ring s, Eq s) => Show (Multivector s a) where
    show (Multivector m)
        | M.size m == 0 = "0"
    show (Multivector m) = intercalate " + " $ fmap showTerm $ M.assocs m
        where
            showTerm (Blade vec, s)
                | s == one = if (length vec == 0) then "1" else show (Blade vec)
                | otherwise = (show s) ++ (show (Blade vec))

newtype Blade a = Blade { vecs :: [a] } deriving (Eq, Ord)

instance Show a => Show (Blade a) where
    show (Blade vec) = concatMap show vec

instance Functor Blade where
    fmap f (Blade xs) = Blade (fmap f xs)

--mostly for testing
--the "mother algebra"
--pick out some of the vecotrs from it to work with, have to be careful about the selection
data Mother = Pos String Int
            | Neg String Int
            deriving (Eq)

instance Ord Mother where
    (Pos _ i) <= (Pos _ j) = i <= j
    (Pos _ _) <= (Neg _ _) = True
    (Neg _ _) <= (Pos _ _) = False
    (Neg _ i) <= (Neg _ j) = i <= j

instance Show Mother where
    show (Pos a _) = a
    show (Neg a _) = a

--sorta dorky that they are the same definition
instance InnerProduct Mother Double where
    dot (Pos _ _) (Neg _ _) = zero
    dot (Neg _ _) (Pos _ _) = zero
    dot (Neg _ i) (Neg _ j)
        | i == j = (neg one)
        | otherwise = zero
    dot (Pos _ i) (Pos _ j)
        | i == j = one
        | otherwise = zero

instance InnerProduct Mother Sym where
    dot (Pos _ _) (Neg _ _) = zero
    dot (Neg _ _) (Pos _ _) = zero
    dot (Neg _ i) (Neg _ j)
        | i == j = (neg one)
        | otherwise = zero
    dot (Pos _ i) (Pos _ j)
        | i == j = one
        | otherwise = zero


--inner product space
class Scalar s where

instance Scalar Double where

instance Scalar Sym where

class Scalar s => InnerProduct a s where
    dot :: a -> a -> s

--just want an abstraction for 0 and 1 numerically, I believe this is a Ring
--i guess the way it will be used means it can be any scalar field for the coefficients
class Ring a where
    zero :: a
    one :: a
    add :: a -> a -> a
    sub :: a -> a -> a
    sub a b = add a (neg b)
    neg :: a -> a
    neg a = sub zero a
    mul :: a -> a -> a
    {-# MINIMAL zero, one, add, mul, (sub | neg)#-}

instance Ring Double where
    zero = 0
    one = 1
    add a b = a + b
    mul a b = a * b
    neg a = negate a

data Sym = Var String
         -- | Zero
         -- | One
         | Lit Double
         | Add Sym Sym
         | Mul Sym Sym
         | Negative Sym
         deriving (Eq)

instance Show Sym where
    show (Var s) = s
    show (Lit d) = show d
    --show Zero = "0"
    --show One = "1"
    show (Add a b) = "(" ++ show a ++ " + " ++ show b ++ ")"
    show (Mul a b) = "(" ++ show a ++ " * " ++ show b ++ ")"
    show (Negative a) = "(" ++ "-" ++ show a ++ ")"

instance Ring Sym where
    zero = Lit 0
    one = Lit 1

    add (Lit 0) b = b
    add a (Lit 0) = a
    add (Lit a) (Lit b) = Lit $ a+b
    add a b = Add a b

    mul (Lit 0) b = (Lit 0)
    mul a (Lit 0) = (Lit 0)
    mul (Lit 1) b = b
    mul a (Lit 1) = a
    mul (Lit a) (Lit b) = Lit $ a*b
    mul (Negative a) (Negative b) = mul a b
    mul (Negative a) b = Negative (mul a b)
    mul a (Negative b) = Negative (mul a b)
    mul a b
        | a == b = mul (Lit 2) a
        | otherwise = Mul a b

    neg (Lit 0) = (Lit 0)
    neg (Negative a) = a
    neg (Lit a) = (Lit (-a))
    neg a = (Negative a)


makeMulti :: (Eq s, Ring s, InnerProduct a s, Ord a) => [(Blade a,s)] -> Multivector s a
makeMulti = Multivector . M.filter (/=zero) . M.fromListWith add

canon :: (Eq s, Ring s, InnerProduct a s, Ord a) => Multivector s a -> Multivector s a
canon (Multivector m) = makeMulti $ fmap canonBlade (M.assocs m)
    where
        canonBlade (Blade vec, c) = case vec of
                                        [] -> (Blade vec, c)
                                        [x] -> (Blade vec, c)
                                        [x,y] -> case compare x y of
                                                    EQ -> (Blade [], mul c (dot x y))
                                                    LT -> (Blade vec, c)
                                                    GT -> (Blade [y,x], neg c)
                                        --due to having the orthonormal assumption this can only ever return a single blade
                                        _ -> let (Blade vec', c') = helper [] vec in (Blade vec', mul c' c)
        -- should do this more in a back to front manner
        helper front [] = (Blade front, one)
        helper front [x] = (Blade $ front ++ [x], one)
        helper front (a:b:xs)
            | a < b = helper (front ++ [a]) (b:xs)
            | otherwise = let (Blade vec', c') = canonBlade (Blade [a,b], one)
                          in (\(a, c) -> (a, mul c c')) $ helper [] (front ++ vec' ++ xs) --have to rescale at the end

--not correct
isBlade :: Multivector s a -> Bool
isBlade (Multivector m) = (1==) $ M.size m

isGrade :: Int -> Multivector s a -> Bool
isGrade n (Multivector m) = let lens = n : fmap (\(Blade vec, _) -> length vec) (M.assocs m)
                            in and (zipWith (==) lens (tail lens))

--finds the grade of a Multivector, can fail
grade :: Multivector s a -> Int
grade (Multivector a) = let lens = fmap (\(Blade vec, _) -> length vec) (M.assocs a)
                        in if and (zipWith (==) lens (tail lens)) then head lens else error "Not a single grade"

extractGrade :: Int -> Multivector s a -> Multivector s a
extractGrade n (Multivector m) = Multivector $ M.filterWithKey (\(Blade vec) _ -> length vec == n) m

unMulti :: Multivector s a -> [(Blade a, s)]
unMulti (Multivector a) = M.assocs a
toMulti :: (Ord a, InnerProduct a s, Ring s, Eq s) => [(Blade a,s)] -> Multivector s a
toMulti = Multivector . M.filter (/=zero) . M.fromListWith add

--are the operations a ring?
mAdd :: (Ring s, InnerProduct a s, Ord a, Eq s) => Multivector s a -> Multivector s a -> Multivector s a
mAdd (Multivector a) (Multivector b) = canon $ toMulti $ (M.assocs a) ++ (M.assocs b)

mSub a b = mAdd a (scale (neg one) b)

mMul :: (Ring s, InnerProduct a s, Ord a, Eq s) => Multivector s a -> Multivector s a -> Multivector s a
mMul (Multivector a) (Multivector b) = canon $ toMulti $ mult <$> (M.assocs a) <*> (M.assocs b)
    where
        mult (Blade vec1, c1) (Blade vec2, c2) = (Blade $ vec1 ++ vec2, mul c1 c2)

outer :: (Ring s, InnerProduct a s, Ord a, Eq s) => Multivector s a -> Multivector s a -> Multivector s a
outer (Multivector a) (Multivector b) = canon $ toMulti $ concat $ outerBlade <$> (M.assocs a) <*> (M.assocs b)
    where
        outerBlade x@(Blade vec1, c1) y@(Blade vec2, c2) = unMulti $ extractGrade ((length vec1) + (length vec2)) (mMul (toMulti [x]) (toMulti [y]))

inner :: (Ring s, InnerProduct a s, Ord a, Eq s) => Multivector s a -> Multivector s a -> Multivector s a
inner (Multivector a) (Multivector b) = canon $ toMulti $ concat $ innerBlade <$> (M.assocs a) <*> (M.assocs b)
    where
        innerBlade x@(Blade vec1, c1) y@(Blade vec2, c2) = unMulti $ extractGrade (abs $ (length vec1) - (length vec2)) (mMul (toMulti [x]) (toMulti [y]))

leftContract :: (Ring s, InnerProduct a s, Ord a, Eq s) => Multivector s a -> Multivector s a -> Multivector s a
leftContract (Multivector a) (Multivector b) = canon $ toMulti $ concat $ lBlade <$> (M.assocs a) <*> (M.assocs b)
    where
        lBlade x@(Blade vec1, c1) y@(Blade vec2, c2) = unMulti $ extractGrade ((length vec2) - (length vec1)) (mMul (toMulti [x]) (toMulti [y]))

rightContract :: (Ring s, InnerProduct a s, Ord a, Eq s) => Multivector s a -> Multivector s a -> Multivector s a
rightContract (Multivector a) (Multivector b) = canon $ toMulti $ concat $ rBlade <$> (M.assocs a) <*> (M.assocs b)
    where
        rBlade x@(Blade vec1, c1) y@(Blade vec2, c2) = unMulti $ extractGrade ((length vec1) - (length vec2)) (mMul (toMulti [x]) (toMulti [y]))

scale s (Multivector m) = canon $ Multivector $ M.map (mul s) m

rev (Multivector a) = canon $ Multivector $ M.mapKeys (\(Blade v) -> (Blade $ reverse v)) a

--haven't worked out a great way of doing an inverse of a multivector (can sometimes fail)
--but for a blade it is straightforward
invBlade a = scale (1/(norm2 a)) (rev a)

norm2 a = let a_dag = rev a
              s = extractGrade 0 (a_dag `mMul` a)
          in case (\(Multivector a) -> M.elems a) s of
                [] -> zero
                [x] -> x
                _ -> error "Uhh... something weird in the squared norm"

--project a onto b
-- works for a : multi, b : blade
project a b = (leftContract a b ) `mMul` (invBlade b)

--reflect a by b
--it is the reflection of all terms in a by b
reflect a b = fromTerms $ fmap (\x -> reflectBlade x b) (toTerms a)
    where
        reflectBlade a' b = let j = grade a'
                                k = grade b
                            in scale ((-1)^(j*(k+1))) (b `mMul` a' `mMul` (invBlade b))

--rotate a by theta in i
rotate a theta i = let z = (toMulti [(Blade [], cos (theta / 2))]) `mAdd` (scale (sin (theta / 2)) i)
                   in (rev z) `mMul` a `mMul` z
