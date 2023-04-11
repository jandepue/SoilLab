import numpy

## READ FILES

def ReadKSAT_0(filename):
    fID=open(filename)
    s=fID.read().replace('\r\n','\n').split('\n')
    nL=len(s)
    
    KsID = 0
    Tnorm = 0
    Ks = 0
    SampleArea = 0
    SampleLength = 0
    
    for iL in range(nL):
        if 'File name'  in s[iL]:
            readID = s[iL].replace('"','').replace(';','').split('\t')[1]
            KsID = readID
        elif ('Ks Soil normalized at' in s[iL]) & ('m/s' in s[iL]) :
#            S=s[iL].replace('Ks Soil normalized at 10,0 \xc2\xb0C [m/s]\t','').replace(',','.').replace(';','').replace('"','')
            Tnorm = float(s[iL].replace('"','').replace('Ks Soil normalized at ','').split(' ')[0].replace(',','.'))
            S=s[iL].split('\t')[1].replace(',','.').replace(';','').replace('"','')
            if 'NaN' in S:
                Ks = numpy.nan
            else:
                Ks = float(S)
        elif 'A_sample  [' in s[iL]:
            SampleArea = float(s[iL].replace('"','').replace(',','.').split('\t')[1])
        elif 'L_sample [cm]' in s[iL]:
            SampleLength = float(s[iL].replace('"','').replace(',','.').split('\t')[1])
    
    # Check
    if  (KsID == 0) |\
        (Ks == 0) |\
        (SampleArea == 0) |\
        (SampleLength == 0):
        raise NameError('Error whilst reading KSAT File')
    else:
        return KsID, Ks, Tnorm, SampleArea, SampleLength

## TOOLS

def CapitalizeArray(X):
    shape = X.shape
    X = numpy.reshape(X,-1)
    nX = X.size
    for iX in range(nX):
        X[iX] = X[iX].upper()
    X = numpy.reshape(X,shape)
    return X

def RemoveInitialZeros(X):
    while (X[0] == '0') & (len(X)>1):
        X = X[1:]
    return X


## KSAT FUNK

def RemoveOutlier(X):
    maxvalue=1000.0/(100*60*60*24) # 1000 cm /d
    n=X.size
    if n > 0:
        X_filtered=X[X<maxvalue]
        n_removed=n-X_filtered.size
        if n_removed > 0:
            print('%s / %s outliers removed'%(n_removed,n))
        if n_removed==n:
            print('\tALERT ALERT: the data is completely outliered, replaced by default value')
#            X_filtered=numpy.array([maxvalue,])
            X_filtered=numpy.array([numpy.nan,])
    else:
#        X_filtered=X
        X_filtered=numpy.array([numpy.nan,])
    return X_filtered

def gstd(X):
    import numpy
    from scipy import stats
#    std=10**(numpy.sqrt(numpy.sum(numpy.log10(X/stats.gmean(X))**2)/numpy.size(X)))
    std=numpy.exp(numpy.sqrt(numpy.sum(numpy.log(X/stats.gmean(X))**2)/numpy.size(X)))
#    std=numpy.exp(numpy.std(numpy.log(X),ddof=1))
    return std

def Sigma_Finney(x,n):
    n=float(n)
    Sigma=1.0
    rel_add=1.0
    i=1
    while rel_add>=1e-2:
        add=(x**i * (n-1.0)**(2*(i-1.0)+1.0))/(n**(i) * numpy.prod(numpy.arange(n-1,n-1+2*i,2))/(n-1) * numpy.math.factorial(i))
        Sigma=Sigma+add
#        print add
        rel_add=add/Sigma
        i=i+1
    return Sigma
    
def gmean_Finney(X,n):
    import numpy
    #    import Sigma_Finney
    varlog=numpy.var(numpy.log(X),ddof=1)
    meanlog=numpy.mean(numpy.log(X))
    mean=numpy.exp(meanlog)*Sigma_Finney(varlog/2,n)
    return mean
    
def gstd_Finney(X,n):
    import numpy
    #    import Sigma_Finney 
    varlog=numpy.var(numpy.log(X),ddof=1)
    meanlog=numpy.mean(numpy.log(X))
    if n>1:
        std=numpy.sqrt(numpy.exp(2*meanlog)*(Sigma_Finney(2*varlog,n)-Sigma_Finney((n-2.0)/(n-1.0)*varlog,n)))
    else:
        std=0.0
    return std
    
    
def gmean_Warrick(X,n):
    import numpy
    #    import Sigma_Finney
    varlog=numpy.var(numpy.log(X),ddof=1)
    meanlog=numpy.mean(numpy.log(X))
    mean=numpy.exp(meanlog+0.5*varlog)
    return mean
    
def gstd_Warrick(X,n):
    import numpy
    #    import Sigma_Finney 
    varlog=numpy.var(numpy.log(X),ddof=1)
    meanlog=numpy.mean(numpy.log(X))
    std=numpy.sqrt(numpy.exp(2*meanlog+varlog)*(numpy.exp(varlog)+1.0))
    return std


