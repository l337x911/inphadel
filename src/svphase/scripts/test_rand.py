import timeit
import numpy as na

BUFFER_SIZE = 10000
N = 50000

def rand():
  return na.random.rand(1)

def rand_gen():
  while True:
    for i in na.random.rand(BUFFER_SIZE):
      yield i

def repeat_rand():
  x = 0
  for i in xrange(N):
    x += rand()

def repeat_rand_gen():
  x = 0
  gen = rand_gen()
  for i in xrange(N):
    x += gen.next()
   
  

if __name__ == '__main__':
  import timeit
  print "RAND"
  print(timeit.repeat("repeat_rand()", setup="from __main__ import repeat_rand", repeat=3, number=50))
  print "RAND GEN"
  print(timeit.repeat("repeat_rand_gen()", setup="from __main__ import repeat_rand_gen", repeat=3, number=50))

