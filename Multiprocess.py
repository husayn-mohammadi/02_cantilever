import multiprocessing as mp
import os
import sys
import openseespy.opensees as ops
import eqsig
import numpy as np
import matplotlib.pyplot as plt
import time

def idaFunction():
    exec(open("IDA.py").read())

# idaFunction()

if __name__ == '__main__':
    
    p1 = mp.Process(target=idaFunction)
    p1.start()
    print("p1 finished")




# result_1 = []
# result_2 = []
# result_3 = []

# def calculation1(numbers):
#     for number in numbers:
#         result_1.append((number ** 3) **0.5)
        
# def calculation2(numbers):
#     for number in numbers:
#         result_2.append((number ** 5) **0.5)
        
# def calculation3(numbers):
#     for number in numbers:
#         result_3.append((number ** 7) **0.5)
        

# if __name__ == '__main__': # when multiprocessing this if command is a MUST
#     number_list = list(range(int(1e6)))
    
    
#     p1 = mp.Process(target = calculation1, args = (number_list, ))
#     p2 = mp.Process(target = calculation2, args = (number_list, ))
#     p3 = mp.Process(target = calculation3, args = (number_list, ))
    
    
#     # Multiprocessing
#     t1 = time.time()
#     p1.start()
#     p2.start()
#     p3.start()
#     t2 = time.time()
#     durMulti = t2-t1
#     print(f"Multiprocessing took {t2-t1} sec")   
    
#     temp1 = result_1
#     temp2 = result_2
#     temp3 = result_3
    
#     # Single Processing
#     t1 = time.time()
#     calculation1(number_list)
#     calculation2(number_list)
#     calculation3(number_list)
#     t2 = time.time()
#     durSingle = t2-t1
#     print(f"Single Processing took {t2-t1} sec")

#     print(f"Multiprocessing was {durSingle/durMulti:.1f}x faster!")

#     print(temp1 == result_1)
#     print(temp2 == result_2)
#     print(temp3 == result_3)





























