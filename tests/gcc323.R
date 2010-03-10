require(spatstat)
data(lansing)

# critical R values that provoke GCC bug #323

a <- marktable(lansing, R=0.25)
a <- marktable(lansing, R=0.21)
a <- marktable(lansing, R=0.20)
a <- marktable(lansing, R=0.10)
