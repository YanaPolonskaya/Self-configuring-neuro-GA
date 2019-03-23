from multiprocessing import  Pool

import sys
{newpath}

def NN(i):

    if i==0:
#GP_Keras_generated_part
        {importer}
        

    #return MyRNN(i)
#End_GP_Keras_generated_part

def Calculate_fitnesses(size_pop):
    #print(sys.path)
    pool=Pool(processes=45)
    #start_time = time.time()
    results = [pool.map(NN,range(size_pop))]
    #print(results[0])
    print("max:"+str(max(results[0])))
    print("ave:"+str(sum(results[0])/size_pop))
    print("min:"+str(min(results[0])))
    pool.close()
    pool.join()
    return results[0]
    # exiting the 'with'-block has stopped the pool

if __name__ == '__main__':
    print("It's main function")
    Calculate_fitnesses("Called from main function")
    # exiting the 'with'-block has stopped the pool
