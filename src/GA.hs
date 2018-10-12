{-# LANGUAGE TypeSynonymInstances, FlexibleInstances #-}
module GA where

import qualified Data.Map.Strict as M

{-
Want to test some stuff with GA and multivector arithmetic
all done assuming that the basis vectors are orthonormal

if the basis was to have arbitrary signature then a "multiplication table" would have to be carried around
    including it would be useful in defining null vectors and other mixed signature spaces
-}

--a vector needs a representation and a way to order it
--the ordering is so that there is a canonical form to always fall back on
data Vector = Vector String Int deriving (Eq)

instance Show Vector where
    show (Vector name _) = name
instance Ord Vector where
    (Vector _ i) <= (Vector _ j) = i <= j
    compare (Vector _ i) (Vector _ j) = compare i j

type Multivector = M.Map [Vector] Double
type InnerProd = M.Map (Vector, Vector) Double

--define an orthonormal frame with dimension n, make G(n,0)
orthonormalPos :: Int -> ([Multivector], InnerProd)
orthonormalPos n = let vec n = Vector ('e':(show n)) n
                       inner = M.fromList [((vec i, vec j), if i==j then 1 else 0) | i <- [1..n], j <- [i..n]]
                       multis = fmap (\x -> makeMulti [([vec x],1)]) [1..n]
                   in (multis,inner)

--orthonormal wher ethe magnitudes are all negative
orthonormalNeg :: Int -> ([Multivector], InnerProd)
orthonormalNeg n = let vec n = Vector ('e':(show n)) n
                       inner = M.fromList [((vec i, vec j), if i==j then (-1) else 0) | i <- [1..n], j <- [i..n]]
                       multis = fmap (\x -> makeMulti [([vec x],1)]) [1..n]
                   in (multis,inner)

--create a mixed signature orthonormal basis
--make it so the positive elements are ordered first
fromSignature :: (Int, Int) -> ([Multivector], InnerProd)
fromSignature (positive, negative) = let vecP n = Vector ('p':(show n)) n
                                         vecM n = Vector ('m':(show n)) (n+positive)
                                         multis = (fmap (\x -> makeMulti [([vecP x],1)]) [1..positive]) ++ (fmap (\x -> makeMulti [([vecM x],1)])) [1..negative]
                                         posVecs = fmap vecP [1..positive]
                                         negVecs = fmap vecM [1..negative]
                                         innerPos = M.fromList [((vecP i, vecP j), if i==j then 1 else 0) | i <- [1..positive], j <- [i..positive]]
                                         innerNeg = M.fromList [((vecM i, vecM j), if i==j then (-1) else 0) | i <- [1..negative], j <- [i..negative]]
                                         innerBetween = M.fromList [((vecP i, vecM j), 0) | i <- [1..positive], j <- [1..negative]]
                                     in (multis, foldr M.union M.empty [innerPos, innerNeg, innerBetween])


one :: Multivector
one = makeMulti [([],1)]
zero :: Multivector
zero = makeMulti []

--definitions for G(2,0)
e1 = Vector "e1" 1
e2 = Vector "e2" 2
innerProdG2 :: InnerProd
innerProdG2 = M.fromList [((e1,e1),1), ((e2,e2),1), ((e1,e2),0)]

--definitions for G(1,1)
ep = Vector "e+" 1
em = Vector "e-" 2
innerProdG11 :: InnerProd
innerProdG11 = M.fromList [((ep,ep),1), ((em,em),-1), ((ep,em),0)]

--definitions for a 2d null space, G(0,0,2)?
n0 = Vector "n0" 1
n = Vector "n_inf" 2
innerProdG002 :: InnerProd
innerProdG002 = M.fromList [((n0,n0),0),((n,n),0),((n0,n),-1)]

multV :: InnerProd -> Vector -> Vector -> Multivector
multV innerProds x y
    | x == y = M.fromList [([], innerProds M.! (x, y))]
    | x <= y = M.fromList [([x,y], 1)]
    | x >= y = makeMulti [([],2*(innerProds M.! (y,x))), ([y,x],-1)]

makeMulti x = M.filter (/=0) $ M.fromListWith (+) x

mScale :: Double -> Multivector -> Multivector
mScale s = M.map (s*)

isBlade :: Multivector -> Bool
isBlade = (1==) . M.size

--check if all the terms have a vector grade equal to n
isGrade :: Int -> Multivector -> Bool
isGrade n m = let lens = n : (fmap (\(vec,_) -> length vec) $ M.assocs m)
              in and $ zipWith (==) lens (tail lens)


--map the canonincalizing function over all the terms seperately
--currently very innefficient, but it works
canon :: InnerProd -> Multivector -> Multivector
canon inner m = M.filter (/=0) $ foldr (M.unionWith (+)) M.empty $ fmap (canon' inner) $ fmap (\x -> makeMulti [x]) (M.assocs m)
    where
        canon' inner m
            | isGrade 0 m = m
            | isGrade 1 m = m
            | isGrade 2 m = let [([x,y],k)] = M.assocs m --have forced it to be one term, so this should be fine
                            in mScale k (multV inner x y)
            | otherwise = let [(vecs, k)] = M.assocs m
                          in mScale k (helper inner [] vecs) --want to factor out coefficient for simplicity
        --the helper scans forward in the vectors and sees if all the terms are in canonical order
        --if all the adjacent pairs of vectors are in canonical form the whole thing is canonical
        --if a pair is not canonical, then get the canonical form, adjust the multivector apprpriately and get the canonical form of the whole thing

        --simplest case have gotten to the end of the vector, so it is in canonical form
        helper inner front [] = makeMulti [(front,1)]
        helper inner front [x] = makeMulti [(front ++ [x],1)] --a single element can't be out of order
        helper inner front (a:b:xs) = let cForm = multV inner a b
                                      in if cForm == makeMulti [([a,b],1)] then --the vectors are in canonical form, move to next pair
                                            helper inner (front ++ [a]) (b:xs)
                                         else --not in canonical form, have to reform the appropriate multivector and then get the canonical form of that
                                            canon inner $ makeMulti $ fmap (\(vec,coe) -> (front ++ vec ++ xs, coe)) (M.assocs cForm)

mAdd :: Multivector -> Multivector -> Multivector
mAdd x y = M.filter (/=0) $ M.unionWith (+) x y
mSub :: Multivector -> Multivector -> Multivector
mSub x y = mAdd x (mScale (-1) y)

mMult :: InnerProd -> Multivector -> Multivector -> Multivector
mMult inner x y = makeMulti $ concat $ ((mult' inner) <$> M.assocs x <*> M.assocs y)
    where
        mult' inner (vec1, coe1) (vec2, coe2) = M.assocs $ canon inner $ makeMulti [(vec1++vec2, coe1*coe2)]

--some more definitions

--I suppose these are only correct for vectors and not general multivectors
--want to define the general equations here
mOuterV :: InnerProd -> Multivector -> Multivector -> Multivector
mOuterV inner x y = mScale (1/2) ((mMult inner x y) `mSub` (mMult inner y x))
mInnerV :: InnerProd -> Multivector -> Multivector -> Multivector
mInnerV inner x y = mScale (1/2) ((mMult inner x y) `mAdd` (mMult inner y x))

{-
The outer and inner products can be defined for blades as
outer a b = extractGrade (i+j) (mult a b) where i = grade a, j = grade b
inner a b = extractGrade (i-j) (mult a b)

so the inner and outer products are
    outerBlabe <$> assocs a <*> assocs b
    innerBlade <$> assocs a <*> assocs b

However it is dubious whether or not to make the inner product extract
    i-j
    or
    |i-j|
-}

mOuter :: InnerProd -> Multivector -> Multivector -> Multivector
mOuter inner x y = makeMulti $ concat $ (outerBlade inner) <$> M.assocs x <*> M.assocs y
    where
        outerBlade inner a@(vecA, coeA) b@(vecB, coeB) = let i = length vecA
                                                             j = length vecB
                                                         in M.assocs $ extractGrade (i+j) (mMult inner (makeMulti [a]) (makeMulti [b]))

mInner :: InnerProd -> Multivector -> Multivector -> Multivector
mInner inner x y = makeMulti $ concat $ (innerBlade inner) <$> M.assocs x <*> M.assocs y
    where
        innerBlade inner a@(vecA, coeA) b@(vecB, coeB) = let i = length vecA
                                                             j = length vecB
                                                         in M.assocs $ extractGrade (abs $ i-j) (mMult inner (makeMulti [a]) (makeMulti [b]))
--extract a certain grade from a multivector
extractGrade :: Int -> Multivector -> Multivector
extractGrade n m = M.filterWithKey (\vec _ -> length vec == n) m

mRev :: InnerProd -> Multivector -> Multivector
mRev inner x = canon inner $ M.mapKeys reverse x

--start defining standard linear mappings
--rotation
--projection
--mirroring


{- good and functioning
data Vector = Vector Int deriving (Eq, Ord)

e1 = basis 1
e2 = basis 2
e3 = basis 3

basis i = makeMulti [([Vector i], 1)]

instance Show Vector where
    show (Vector i) = "e"++ (show i)

--want to sort with and count number of swaps
sortWithSwaps xs = gnome 0 [] xs
    where
        gnome s front [] = (front, s)
        gnome s front [a] = (front ++ [a], s)
        gnome s [] (a:b:xs)
            | a <= b = gnome s [a] (b:xs)
            | otherwise = gnome (s+1) [b] (a:xs)
        gnome s front (a:b:xs)
            | a <= b = gnome s (front ++ [a]) (b:xs)
            | otherwise = gnome (s+1) (init front) (last front : b : a : xs)

--a multivector stores the blades ( the base vector and the coefficient )
type Multivector = M.Map [Vector] Double

--the canonical form of a blade has the vectors sorted, negates the coefficient as needed, and cancels the needed terms (e1e1 = 1)
--sort of have these backwards and it makes a few spots jumbly later, likely want canon (vecs,coe) -> (vecs, coe)
canon x@([],coe) = x
canon (vecs, coe) = let (vecs', s) = sortWithSwaps vecs
                        vecs'' = foldr (\x acc -> if null acc then [x] else if x == (head acc) then tail acc else x:acc) [] vecs'
                    in if even s then ( vecs'', coe) else ( vecs'', -coe)

--make a multivector and make sure it is in canonical form
makeMulti :: [([Vector], Double)] -> Multivector
makeMulti xs = M.filter (/=0) $ M.fromListWith (+) (fmap canon xs)

--adding multivectors is then pretty easy as long as all the blades are in the canonical form
mAdd :: Multivector -> Multivector -> Multivector
mAdd x = M.filter (/=0) . M.unionWith (+) x

mSub x = M.filter (/=0) . M.unionWith (-) x

--multiplication is then multiplication between all of the individual terms summed together
mMul :: Multivector -> Multivector -> Multivector
mMul a b = makeMulti $ mult <$> (M.assocs a) <*> (M.assocs b)
    where
        mult (b1, c1) (b2, c2) = canon ((b1++b2),(c1*c2))

--the reverse of a multivector, reverse all the basis of the blades and the nput it back in canonical form
mRev :: Multivector -> Multivector
mRev = makeMulti . M.assocs . M.mapKeys reverse

mScale :: Double -> Multivector -> Multivector
mScale s = M.map (s*) -- == mMul (makeMulti [([],s)])

--the inner product is defined by a.b = 1/2*(ab+ba)
mInner :: Multivector -> Multivector -> Multivector
mInner a b = mScale (1/2) ((a `mMul` b) `mAdd` (b `mMul` a))

--the outer product is a^b = 1/2*(ab-ba)
mOuter :: Multivector -> Multivector -> Multivector
mOuter a b = mScale (1/2) ((a `mMul` b) `mSub` (b `mMul` a))

--some extra and useful manipulations and inspections
extractGrade :: Int -> Multivector -> Multivector --extract the given grade elements
extractGrade n m = M.filterWithKey (\k _ -> length k == n) m

isBlade :: Multivector -> Bool --whether or not all the terms are the same grade or not
isBlade m = let vecs = fmap (length . fst) $ M.assocs m in and (zipWith (==) vecs (tail vecs))

--isGrade n m -> extractGrade n m == m
isGrade n m
    | not $ isBlade m = False --can it have a grade if it isn't a blade?
    | otherwise = (length $ fst $ head $ M.assocs m) == n


--ehh?
instance Num Multivector where
    (+) = mAdd
    (-) = mSub
    (*) = mMul
    negate = mScale (-1)
    abs x = M.map sqrt $ x * (mRev x)
    signum = undefined -- ?
    fromInteger i = makeMulti [([],fromIntegral i)]

-}








{- Old not great

data Base = Scalar Double | Vector Int deriving (Eq)

instance Show Base where
    show (Scalar d) = show d
    show (Vector i) = "e"++show i

data Vector = E Int deriving (Eq)
instance Show Vector where
    show (E i) = "e"++(show i)

data Expr a = Inner a a
            | Outer a a
            | Geo a a
            | Sum a a
            | Val a
            deriving (Eq)

instance Show a => Show (Expr a) where
    show (Inner a b) = "(" ++ show a ++ "." ++ show b ++ ")"
    show (Outer a b) = "(" ++ show a ++ "^" ++ show b ++ ")"
    show (Geo a b) = "(" ++ show a ++ show b ++ ")"
    show (Sum a b) = "(" ++ show a ++ "+" ++ show b ++ ")"
    show (Val a) = show a

eval :: Expr Base -> Expr Base
eval (Inner (Vector i) (Vector j))
    | i == j = Val $ Scalar 1
    | i /= j = Val $ Scalar 0
eval (Outer (Vector i) (Vector j))
    | i == j = Val $ Scalar 0
eval (Geo (Scalar a) (Scalar b)) = Val $ Scalar (a*b)
eval (Sum (Scalar a) (Scalar b)) = Val $ Scalar (a+b)
eval (Sum (Scalar 0) a) = eval $ Val a
eval (Sum a (Scalar 0)) = eval $ Val a
eval (Geo (Scalar 0) a) = Val $ Scalar 0
eval (Geo a (Scalar 0)) = Val $ Scalar 0
eval (Geo (Scalar 1) a) = eval $ Val a
eval (Geo a (Scalar 1)) = eval $ Val a
eval x = x

e1 = Vector 1
e2 = Vector 2
-}
