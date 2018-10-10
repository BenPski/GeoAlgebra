# GA

This is for messing around with geometric algebra calculations.

Works by first creating the generating vectors and defining the inner product relations, then all the computations can be determined from those by picking a canonical ordering of the vectors.

To work with G(3,0) can do:
```haskell
([e1,e2,e3],innerProd) = orthonormal 3
```
That will define the the 3 basis vectors and the inner product mappings. Then, use *mAdd* and *mMult* to do addition and the geometric product.

For something like G(1,1), a 2D mixed signature:
```haskell
--generating vector definition
ep_base = Vector "e+" 1
em_base = Vector "e-" 2

--the basis vectors
ep = makeMulti [([ep_base],1)]
em = makeMulti [([em_base],1)]

--the inner products
inner = Map.fromList [((ep_base, ep_base),1), ((em_base, em_base),-1), ((ep_base, em_base),0)]
```

Or a 2D null space, G(0,0,2)?
```haskell
--generating vector definition
e0_base = Vector "e0" 1
einf_base = Vector "e" 2

--the basis vectors
e0 = makeMulti [([e0_base],1)]
e = makeMulti [([einf_base],1)]

--the inner products
inner = Map.fromList [((e0_base, e0_base),0), ((einf_base, einf_base),0), ((e0_base, einf_base),-1)]
```

Then things like addition don't need the inner product definition, but multiplication does.
```haskell
mAdd x y
mMult inner x y
```
