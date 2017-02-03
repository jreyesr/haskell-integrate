import System.IO
import System.Environment

prompt :: String -> IO String
prompt text = do
    putStr text
    hFlush stdout
    getLine

approx :: Fractional a => (a -> a) -> [a] -> [a] -> a
approx f xs ws = sum [w * f x | (w,x) <- zip ws xs]

integrateNoEnds :: Fractional a => a -> [a] -> (a -> a) -> a -> a -> Int -> a
--Sum the contributions of every xs[i], multiplied by its weight ws[i]
integrateNoEnds v vs f a b n = approx f xs ws * h / v where
  --Divide [a,b] in N(vs)*n divisions -> m
  m = fromIntegral (length vs) * n
  --Find the length of a division -> h
  h = (b-a) / fromIntegral m
  --Replicate the weights n times
  ws = concat $ replicate n vs
  --Find the values of x for evaluation
  xs = [a + h/2 + h * fromIntegral i | i <- [0..m-1]]

integrateEnds :: Fractional a => a -> [a] -> (a -> a) -> a -> a -> Int -> a

--Sum the contributions of every xs[i], multiplied by its weight ws[i]
integrateEnds v vs f a b n = approx f xs ws * h / v where
  --Divide [a,b] in (N(vs)-1)*n divisions -> m
  m = fromIntegral (length vs - 1) * n
  --Find the length of a division -> h
  h = (b-a) / fromIntegral m
  --Find the weights -> ws
  ws = overlap n vs
  --Find the values of x for evaluation
  xs = [a + h * fromIntegral i | i <- [0..m]]

--Repeats xs n times, overlapping ONE element
overlap :: Num a => Int -> [a] -> [a]
overlap n [] = []
overlap n (x:xs) = x:inter n xs where
  inter 1 ys     = ys
  inter n []     = x:inter (n-1) xs
  inter n [y]    = (x+y):inter (n-1) xs
  inter n (y:ys) = y:inter n ys

--Definitions of common integrals
--integrateLeft: Includes ends, left one weighs 1, left one weighs 0 (nothing)
integrateLeft = integrateEnds 1 [1,0]
--integrateRight: Same, but on reverse
integrateRight = integrateEnds 1 [0,1]
--integrateMiddle: Ends not included, only one point (middle) weighing 1
integrateMiddle = integrateNoEnds 1 [1]
--integrateTrapezium: Includes ends, weighing the same, must divide by 2
integrateTrapezium = integrateEnds 2 [1,1]
--integrateSimpson: Includes ends AND middle, weighing 1, 1 and 4
integrateSimpson = integrateEnds 3 [1,4,1]

--Definition of adaptive integral (no divisions required)
integrateAdaptive method f a b epsilon =
  --If div10 is "close enough" to div5, use div10, else:
  if abs (div10 - div5) < e then div10 else
    --Divide in two parts, [a,middle] y [middle,b]. Recurse
    integrateAdaptive method f a middle (Just e) + integrateAdaptive method f middle b (Just e)
          --Normally calculated, with 5 divisions -> div5
    where div5 = method f a b 5
          --Same with 10 divisions -> div10
          div10 = method f a b 10
          middle = (a + b) / 2
          e = maybe 1e-5 id epsilon

integrate2 method f a b c d divX divY = method (\y -> integratePartial2 method f a b y divX) c d divY
integratePartial2 method f a b k n = method (\x -> f x k) a b n

integrate3 method f a b c d e g divX divY divZ = method (\z -> integratePartial3 method f a b c d z divX divY) e g divZ
integratePartial3 method f a b c d k divX divY = integrate2 method (\x y -> f x y k) a b c d divX divY

integrateAdaptive2 method f a b c d epsilon = integrateAdaptive method (\y -> integrateAdaptivePartial2 method f a b y epsilon) c d epsilon
integrateAdaptivePartial2 method f a b k epsilon = integrateAdaptive method (\x -> f x k) a b epsilon

integrateAdaptive3 method f a b c d e g epsilon = integrateAdaptive method (\z -> integrateAdaptivePartial3 method f a b c d z epsilon) e g epsilon
integrateAdaptivePartial3 method f a b c d k epsilon = integrateAdaptive2 method (\x y -> f x y k) a b c d epsilon

integrateGaussLegendre f a b n = d*sum [ w x*f(m + d*x) | x <- roots ]
  where d = (b - a)/2
        m = (b + a)/2
        w x = 2/(1-x^2)/(legendreP' n x)^2
        roots = map (findRoot (legendreP n) (legendreP' n) . x0) [1..n]
        x0 i = cos (pi*(i-1/4)/(n+1/2))

legendreP n x = go n 1 x
  where go 0 p2 _  = p2
        go 1 _  p1 = p1
        go n p2 p1 = go (n-1) p1 $ ((2*n-1)*x*p1 - (n-1)*p2)/n

legendreP' n x = n/(x^2-1)*(x*legendreP n x - legendreP (n-1) x)

findRoot f df = fixedPoint (\x -> x - f x / df x)

fixedPoint f x | abs (fx - x) < 1e-15 = x
               | otherwise = fixedPoint f fx
  where fx = f x

targetFunc2 x y = recip $ (1+x^2+y^2)^2
targetFunc x y z = recip $ (1+x^2+y^2+z^2)^2

low = (-10)
high = 10

main = do
  (divsStr: _) <- getArgs
  let divs = read divsStr
--  print $ integrate3 integrateSimpson targetFunc low high low high low high divsX divsY divsZ
  print $ integrate3 integrateGaussLegendre targetFunc low high low high low high divs divs divs
--  print $ integrateMiddle3 targetFunc (-100) 100 (-100) 100 (-100) 100 750 750 750
--  print $ integrateAdaptive3 integrateSimpson targetFunc (-100) 100 (-100) 100 (-100) 100 1e-4
