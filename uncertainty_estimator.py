import time
import argparse
import numpy as np
import pandas as pd
from scipy.special import zeta
from functools import partial
from multiprocessing import Pool, cpu_count

def init():
    global rng
    rng = np.random.default_rng()

def random_walk(d,n_sample):
    """
        Random walk procedure to estimate the extrapolation variable and its uncertainty. The function is set up so it can be used with Pool.map function 
        to run in parallel.
        
        Parameters
        ----------
        d : (tuple) Contains data with: 
                                        ext   : The list of starting extrapolated values
                                        thrs  : Threshold for stopping the procedure
                                        n     : The number of maximal iterations
        n_sample : (int) Number of walks per CPU 
    """
    ext, thrs, n = d #unpack the data
    guesses = []
    rng = np.random.default_rng() # create a new RNG generator to not share the same across multiple processes

    # Prepare pool of random numbers between 0 and 1 to not call random generator every time. 
    # Number prepared samples is dependant of number of walks but reduced by a int((log10(n_sample) -1)  to reduce memory footprint
    samples = rng.uniform(0,1,n_sample//int(np.log10(n_sample)-1.)) 
    _i = 0
    for _ in range(n_sample):
        err = 1
        data = ext[:] #make local copy of extrapolation data
        pm = 1 #bracket 
        while abs(pm) > thrs and len(data) < n:
            pm = abs(data[-1] - data[-2])
            new_val = 2. * pm * samples[_i] + data[-1] - pm 
            data.append(new_val)
            _i += 1
            if _i >= len(samples):
                samples = rng.uniform(0,1,n_sample//int(np.log10(n_sample)-1.))
                _i = 0
        if len(data) == n:
            print("")
            errMsg = f": prematurely ending at iteration #.{len(data)} with current err. {err:12.6e}|{thrs:12.6e} and last extr. value {data[-1]:12.7f}. Increase of N is recommended."
            raise RuntimeError(errMsg)
        guesses.append(data[-1])
    return guesses

class Estimator():
    """
        Estimator class to estimate the uncertainty of the extrapolation through the random walk approach DOI:?????
        
    """
    
    def __init__(self, CPU = 0, N_walk = 1000000, N = 300, thrs = 1e-8, verbose = 0):
        """
            Initialization method

            Parameters
            ----------
            CPU     : (int) Number of CPU cores to use for the procedure. Default is half of all available cores
            N_walk  : (int) Number of random walks to estimated the extrapolation and its uncertainty. Default is 1 000 000.
            N       : (int) Number of iteration for each random walk. Default is 300.
            thrs    : (float) Threshold for stopping the iteration for each random walk. Default is 1e-8.
            verbose : (int) How much to print during procedure. 0 := do not print anything, 1 := print final table, 2:= print used options, datapoints, extrapolated results, final table. 
                            Default is 0.
        """
        
        self.CPU     = CPU 
        if CPU == 0:
            self.CPU = int(cpu_count()*0.5) 
        self.N_walk  = N_walk
        self.N       = N
        self.thrs    = thrs
        self.verbose = verbose
                
    def estimate(self,v,k="", **kwargs):
        """ 
            Function to run the actuall random walks
            
            Parameters:
            -----------
            v      : (tuple) Pair of lists. The first list should contain variables to extrapolated, the second list the corresponding cardinal numbers of basis sets used to estimated variables in the first list.
                             If second list contains 'inf' or 'INF' string at the last place, the corresponding variable is assumed to be an exact value. 
            k      : (str)[optional] The optional header to be printed for 'verbose > 1'. Default is no header
            kwargs : (dict)[optional] Dictionary containing additional info for the procedure. 
                                      'n' : The exponent for the default Truhlar or Riemann extrapolation. Default is 3.
                                      'extr_method': The extrapolation method. Possible options are 'truhlar' to use Truhlar extrapolation, 'riemann' to use Riemann extrapolation or
                                                     function containing user defined extrapolation function. The custom function needs to be in a form of func(Lx,Ex,kwargs). Default is Truhlar extrapolation.
        """
        class Result():
            """
                Result class containg information about the extrapolation results and its estimated uncertainty
                
                Initialization variables
                ------------------------
                Lx    : The array of basis set cardinal numbers used in extrapolations. The array has to be sorted in ascending order.
                Ex    : The array of properties calculated using the basis sets with corresponding cardinal numbers defined in Lx.
                e_x   : The array of extrapolated results. 
                exact : The optional exact value of propertie to compare the extrapolation results with.
                
                Public member variables
                -----------------------
                mean   : The array containing the estimated mean values obtained from random walk procedure.
                std    : The array containing the estimated standard deviations obtained from random walk procedure.
                sigmas : The array of list of estimated sigmas (default size is 3) obtained from random walk procedure.
                
                Private member variable
                -----------------------
                _df : Pandas dataframe for the printing information 
            """
            def __init__(self, Lx, Ex, e_x, exact = None):
                """
                    Initialization method
                    
                    Parameters
                    ----------
                    Lx    : (list)
                    Ex    : (list)
                    e_x   : (list)
                    exact : (float) [optional]
                """
                self.Lx = Lx
                self.Ex = Ex
                self.e_x = e_x
                self.exact = exact

                # assert if we the list lengths agree and we have consistent number of datapoints
                assert len(self.Lx) == len(self.Ex)
                assert len(self.e_x) == len(self.Ex) - 1
                
                if np.array_equal(np.arange(self.Lx[0],self.Lx[-1]), self.Lx):
                    raise ValueError("Data not sorted acccording to their cardinal number")

                self.mean = []
                self.std = []
                self.sigmas = []

                self._df = None
                
            def append(self, mean, std, sigmas):            
                """ 
                    Append method to add results of the random walk 
                    
                    Parameters
                    ----------
                    mean   : estimated mean
                    std    : estimated standard deviation
                    sigmas : estimated sigmas
                """
                self.mean.append(mean)
                self.std.append(std)
                self.sigmas.append(sigmas)
                
            def print(self):        
                """
                    Function to 'nicely' print results using pandas dataframe class
                """
                assert len(self.mean) == len(self.e_x) - 1 
                data_to_print = []
                for i,data in enumerate(zip(self.mean,self.std,self.sigmas)):
                    mu,std,sigma = data
                    sigma1,sigma2,sigma3 = sigma

                    if self.exact is not None:
                        data_to_print.append( [f" {{{self.Lx[i]:2d},{self.Lx[i+1]:2d}}} {{{self.Lx[i+1]:2d},{self.Lx[i+2]:2d}}}",
                                               f"{mu:.7f}",
                                               f"{self.e_x[i+1]:.7f}",
                                               f"{self.exact:.7f}",
                                               f"{std:.7f}",
                                               f"{sigma1:.7f}",
                                               f"{sigma2:.7f}",
                                               f"{sigma3:.7f}",
                                               f"{abs(self.e_x[i+1]-self.Ex[i+2]):.7f}",
                                               f"{abs(self.e_x[i+1]-self.e_x[i]):.7f}",
                                               f"{self.exact-self.e_x[i+1]:.7f}"]
                                            )
                    else:
                        data_to_print.append( [f" {{{self.Lx[i]:2d},{self.Lx[i+1]:2d}}} {{{self.Lx[i+1]:2d},{self.Lx[i+2]:2d}}}",
                                               f"{mu:.7f}",
                                               f"{self.e_x[i+1]:.7f}",
                                               f"{std:.7f}",
                                               f"{sigma1:.7f}",
                                               f"{sigma2:.7f}",
                                               f"{sigma3:.7f}",
                                               f"{abs(self.e_x[i+1]-self.Ex[i+2]):.7f}",
                                               f"{abs(self.e_x[i+1]-self.e_x[i]):.7f}"]
                                            )

                if self.exact is not None:
                    columns = ["Basis sets L","    Est. value","e_x(L)"," Exact"," St.Dev.","1σ","2σ","3σ"," e_x(L)-E(L)"," e_x(L)-e_x(L-1)"," Diff. from exact"]
                else:
                    columns = ["Basis sets L","    Est. value","e_x(L)"," St.Dev.","1σ","2σ","3σ"," e_x(L)-E(L)"," e_x(L)-e_x(L-1)"]    
                d = np.asarray(data_to_print)                
                self._df = pd.DataFrame(d,columns=columns)
                print("-"*150)
                with pd.option_context('display.max_rows', 300, 'display.max_columns', 12, 'display.expand_frame_repr', False): #pandas sometime shorten the printed tables
                    print(self._df)
                print(" ")  

 
        if self.verbose > 1: 
            print(f"CPU               = {self.CPU}")
            print(f"N_walk            = {self.N_walk:,}")
            print(f"thrs              = {self.thrs}")
            print(f"N of extr. steps  = {self.N}")
            print("")
        

        if self.verbose > 1 and len(k) > 0:
            print("Running extrapolation from file:",k)
            
        Ex, Lx = v
        exact = None    
        if Lx[-1] == np.inf or (isinstance(Lx[-1],str) and 'INF' == Lx[-1].upper()):
            exact = Ex[-1]
            Ex = Ex[:-1]
            Lx = Lx[:-1]
            if self.verbose > 1:
                print("Found exact value",exact)        
    
        # There has to be at least 3 datapoints to get two extrapolation results to start the random walk procedure
        assert len(Ex) >= 3 
        assert len(Ex) == len(Lx)

        if self.verbose > 1:
            print("Data to extrapolate")    
            print(Ex)
            print("Basis set cardinal numbers")
            print(Lx)
        
        # default extrapolation exponent is 3
        n_extr = 3
        if kwargs.get("n") is not None:
             n_extr = kwargs.get("n")
                
        # get the extrapolation method
        extr_method = kwargs.get("extr_method") 
        if callable(extr_method):
            ex = extr_method(Lx,Ex,**kwargs)            
            if self.verbose > 1 :
                print("Using custom extrapolation to kickstart",ex)
        elif extr_method is None or extr_method.upper() == "TRUHLAR":
            def truhlar(n0,n1,n):
                l1,e1 = n1
                l1_n = l1**n
                l0,e0 = n0
                l0_n = l0**n
                cbs = e1*l1_n - e0*l0_n
                cbs /= l1_n - l0_n
                return cbs
            ex = [truhlar([Lx[i],Ex[i]],[Lx[i+1],Ex[i+1]],n_extr) for i in range(len(Lx)-1)]
            if self.verbose > 1 :
                print(f'Using Truhlar extrapolation (n={n_extr}) to kickstart',ex)
        elif extr_method.upper() == "RIEMANN":
            ex = [Ex[i] + (Ex[i] - Ex[i-1])*(Lx[i]**(n_extr+1))*(zeta(n_extr+1) - sum([l**(-(n_extr+1)) for l in range(1,Lx[i]+1)])) for i in range(1,len(Lx))]
            if self.verbose > 1 :
                print('Using Riemann extrapolation (n={n_extr}) to kickstart',ex)        
        else:
            raise NotImplementedError("Not implemented method"+str(extr_method))
            
        # prepare the Result object
        results = Result(Lx,Ex,ex,exact)   
        
        # This loops runs over all possible combinations of Lx. 
        # For example, if we have Lx = [2,3,4,5], it will run two iterations:
        # 1st iteration using extrapolations with [2,3], [3,4] 
        # 2nd iteration using extrapolations with [3,4], [4,5]        
        for i in range(0,len(ex)-1):
            t0 = time.time()
            if self.verbose > 1:
                print(f"Random walk from {Lx[i]:2d},{Lx[i+1]:2d}Z and {Lx[i+1]:2d},{Lx[i+2]:2d}Z ",end=" | ")
            pal_random = partial(random_walk, (ex[:i+2], self.thrs, self.N)) # prepare the datapack
            with Pool(self.CPU) as pool:
                guesses = np.asarray(pool.map(pal_random,[self.N_walk//self.CPU]*self.CPU)).flatten()              

            # estimated mean and std of the results
            mu = np.mean(guesses)
            std = np.std(guesses)

            # prepare interpolated percentiles of ~ 1 sigma
            quantiles = np.percentile(guesses,[2.5,97.5])
            
            # subtract the mean value to obtain the errors
            guesses = sorted( abs(guesses - mu) )
            # now estimate the sigmas equal to ~ 68.27%. 95.45% 99.73% from the errors.            
            sigma1 = guesses[int(len(guesses)*0.6827)]
            sigma2 = guesses[int(len(guesses)*0.9545)]
            sigma3 = guesses[int(len(guesses)*0.9973)]
            
            # clean up the memory
            guesses = None
            t1 = time.time() 
            if self.verbose > 1:
                print(f"Mean: {mu:.7f} 2.5% of data lies below {quantiles[0]:.7f} and 97.5% below {quantiles[1]:.7f} | 68.27% {sigma1:.7f} 95.45% {sigma2:.7f} 99.73% {sigma3:.7f} Finished in {t1-t0:.3f}")
            # append this iteration to final results
            results.append(mu,std,[sigma1,sigma2,sigma3])
        if self.verbose:
            results.print()
        return results

    
if __name__ == "__main__":
    
    """ 
        Parameters
        ----------
        i : (str) The filepath to the file containing pairs of cardinal numbers and values to be extrapolated. If 'inf' string in the place of cardinal number, 
                        the value is taken as an exact value.
                        An example of input file is:
                         2    -0.0400183970
                         3    -0.0411736630
                         4    -0.0415978080
                         5    -0.0417856800
                         6    -0.0418812960
                         7    -0.0419349210
                         inf  -0.042044381
        CPU      : (int) Number of CPU cores to use for the procedure. Default is half of all available cores
        N_walk   : (int) Number of random walks to estimated the extrapolation and its uncertainty. Default is 1 000 000.
        N        : (int) Number of iteration for each random walk. Default is 300.
        thrs     : (float) Threshold for stopping the iteration for each random walk. Default is 1e-8.
    """    
    
    default_cpu = cpu_count()//2
    parser = argparse.ArgumentParser(description="Random walk extrapolation")
    parser.add_argument("--i","--input", type=str, required=True, help="Input file to load the datapoints from. The format of file is assumed to contain pairs of cardinal numbers and values on each separate line")
    parser.add_argument("--CPU", type=int, default=default_cpu, help="Number of CPU cores to use (default: half of all available cores).")
    parser.add_argument("--N_walk", type=int, default=1000000, help="Number of walks (default: 1000000).")
    parser.add_argument("--thrs", type=float, default=1e-8, help="Threshold value for random walk (default: 1e-8).")
    parser.add_argument("--N", type=int, default=300, help="Number of extrapolation steps in random walk (default: 300).")

    opts = parser.parse_args()
    
    Lx = []
    Ex = []
    with open(opts.i,'r') as fin:
        for line in fin:                            
            line = line.split()
            if line[0].upper() == "INF":
                Lx.append(line[0])
            else:
                Lx.append(int(line[0]))
            Ex.append(float(line[1]))
    v = [Ex, Lx]          
    estim = Estimator(CPU = opts.CPU, N_walk = opts.N_walk, N = opts.N, thrs = opts.thrs, verbose=2)
    estim.estimate(v,opts.i)

