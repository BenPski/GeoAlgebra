module GA where

import qualified Data.Map.Strict as M
import Data.List (intercalate)

{-

This time around want to avoid having to carry around the inner product space definition
To do this just have to require an inner product space be defined on the base vectors

So the layout goes something like:
    define the generating vectors for the space
    define the inner product space for the generating vectors
    the arithmetic simply follows

also making it a bit clearer that a multivector is a sum of blades with an associated coefficient

issue is now that defining the vector space isn't quite as simple
    no more saying it is orthonormal R^n, at least for now

I think the canonical form stuff doesn't necessarily work for null vectors
but it is fine for the orthonormal ones

-}

--inner product space
class InnerProduct a where
    dot :: a -> a -> Double

--the canonical form of a term
--sorting and trasforming via the inner product
canon :: (InnerProduct a, Ord a) => (Blade a, Double) -> [(Blade a, Double)]
canon a@(Blade vec, coe) = case vec of
                            [] -> [a] --no way it could be out of order
                            [x] -> [a]
                            [x,y] -> blades $ removeZeros $ scale coe (multV x y) --simply do the multiplication
                            --have to step through it
                            _ ->  removeZeros $ scale coe $ helper [] vec
    where
        multV a b
            | a == b = [([], dot a b)]
            | a < b = [([a,b], 1)]
            | a > b = [([],2*(dot b a)), ([b,a], -1)]
        removeZeros = filter (\(_,c) -> c /= 0)
        scale coe = fmap (\(v,c) -> (v, c*coe))
        blades = fmap (\(v,c) -> (Blade v,c))

        helper front [] = [(Blade front, 1)] --completed
        helper front [x] = [(Blade $ front ++ [x], 1)]
        helper front (a:b:xs)
            | a < b = helper (front ++ [a]) (b:xs) --already in order
            | otherwise = let res = multV a b --get the actual form
                              terms = fmap (\(v,c) -> (Blade $ front ++ v ++ xs, c)) res
                          in concatMap canon terms


--something like this
--given a coefficient and the underlying blade, get it into canonical form
--may have to create multiple blades though, like in a null space
    --a nullspace can have ab = 2*(a.b) - ba
    --as the canonical form of the term
--in an orthonormal system that should never happen


--a blade is the free monoid on the generating vectors, but it has to be resricted by the inner product
--so it is isomorphic to a list in its definition, but in applications duplicate items have to be removed in a formulaic way
newtype Blade a = Blade { vecs :: [a] } deriving (Eq, Ord)

instance Show a => Show (Blade a) where
    show (Blade vec) = concatMap show vec

instance Functor Blade where
    fmap f (Blade xs) = Blade (fmap f xs)

--a multivector
newtype MVector a = MVector (M.Map (Blade a) Double) deriving (Eq, Ord)

instance Show a => Show (MVector a) where
    show (MVector m)
        | M.size m == 0 = "0"
    show (MVector m) = intercalate " + " $ fmap showTerm $ M.assocs m
        where
            showTerm (blade, 1)
                | (length $ vecs blade) == 0 = "1"
                | otherwise = (show blade)
            showTerm (blade, s) = (show s) ++ (show blade)

isBlade :: MVector a -> Bool
isBlade (MVector m) = (1==) $ M.size m

isGrade :: Int -> MVector a -> Bool
isGrade n (MVector m) = let lens = n : fmap (\(Blade vec,_) -> length vec) (M.assocs m) in and $ zipWith (==) lens (tail lens)

extractGrade :: Int -> MVector a -> MVector a
extractGrade n (MVector m) = MVector $ M.filterWithKey (\(Blade vec) _ -> length vec == n) m

unMulti :: MVector a -> [(Blade a, Double)]
unMulti (MVector a) = M.assocs a
toMulti :: (Ord a, InnerProduct a) => [(Blade a, Double)] -> MVector a
toMulti = MVector . M.filter (/=0) . M.fromListWith (+)

add :: (Ord a, InnerProduct a) => MVector a -> MVector a -> MVector a
add (MVector a) (MVector b) = toMulti $ (concatMap canon $ M.assocs a) ++ (concatMap canon $ M.assocs b)

sub :: (Ord a, InnerProduct a) => MVector a -> MVector a -> MVector a
sub a b = add a (scale (-1) b)

mult :: (Ord a, InnerProduct a) => MVector a -> MVector a -> MVector a
mult (MVector a) (MVector b) = toMulti $ concat $ mult <$> M.assocs a <*> M.assocs b
    where
        mult (Blade vecs1, coe1) (Blade vecs2, coe2) = canon $ (Blade (vecs1 ++ vecs2), coe1*coe2 )

outer :: (InnerProduct a, Ord a) => MVector a -> MVector a -> MVector a
outer (MVector a) (MVector b) = toMulti $ concat $ outerBlade <$> (M.assocs a) <*> (M.assocs b)
    where
        outerBlade x@(Blade vec1, coe1) y@(Blade vec2, coe2) = unMulti $ extractGrade ((length vec1)+(length vec2)) (mult (toMulti [x]) (toMulti [y]))

inner :: (InnerProduct a, Ord a) => MVector a -> MVector a -> MVector a
inner (MVector a) (MVector b) = toMulti $ concat $ innerBlade <$> (M.assocs a) <*> (M.assocs b)
    where
        innerBlade x@(Blade vec1, coe1) y@(Blade vec2, coe2) = unMulti $ extractGrade ((length vec1)-(length vec2)) (mult (toMulti [x]) (toMulti [y]))

scale s (MVector a) = MVector $ M.map (*s) a

leftContract (MVector a) (MVector b) = toMulti $ concat $ innerBlade <$> (M.assocs a) <*> (M.assocs b)
    where
        innerBlade x@(Blade vec1, coe1) y@(Blade vec2, coe2) = unMulti $ extractGrade ((length vec1)-(length vec2)) (mult (toMulti [x]) (toMulti [y]))
