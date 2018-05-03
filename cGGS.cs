using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Accord.Math;
using Accord.Math.Decompositions;

namespace GGS
{
    class cGGS
    {
        public void fGGS(out List<int> breaks, out List<Double> plotPoints, Double[][] data, int kMax, Double lamb,
            int[] features = null, bool verbose = false)
        {
            //Assign outputs
            breaks = new List<int>();
            plotPoints = new List<double>();//plotPoints was defined for storing likelihoods
            data = cStatistic.fReshape(data);

            List<List<int>> breakPoints = new List<List<int>>();//breakPoints was defined for storing list of breakpoints in each step
            int m, n;//Size of Row and Colomn in data
           
            if(features!=null)
            {
                int nDataLen=data.Length;
                int nFeatureLen=features.Length;
                Double[][] featureData = new Double[nDataLen][];
                for (int ind1 = 0; ind1 < nDataLen; ind1++)
                {
                    featureData[ind1] = new Double[nFeatureLen];
                    for(int ind2=0;ind2<nFeatureLen;ind2++)
                    {
                        featureData[ind1][ind2] = data[ind1][features[ind2]];
                    }
                }
                data = featureData;
            }
           
            m = data.Length;
            n = data[0].Length;
            
            breaks.Add(0);
            breaks.Add(m + 1);

            breakPoints.Add(breaks);

            Double likelihood = fCalculateLikelihood(data, breaks, lamb);
            plotPoints.Add(likelihood);
           
            for (int z = 0; z < kMax; z++)
            {
                int numBreaks = z + 1;
                int ind = -1, newInd = -1;
                Double val = +1, newVal = +1;
                for (  int i = 0; i < numBreaks; i++)
                {
                    List<Double[]> tempData = new List<double[]>();
                    for (int j = breaks[i]; j <= breaks[i + 1] && j<m; j++)
                        tempData.Add(data[j]);
                    fAddBreak(tempData.ToArray(), lamb,out ind,out val);
                    
                    if (val < newVal)
                    {
                        newInd = ind + breaks[i];
                        newVal = val;
                    }//------------------------------if(val<newVal)
                }//------------------------------for(int i=0;i<numBreaks;i++)
                //-----------
                //Check if our algorithm is finished
                if (newVal == 0)
                {
                    return;
                }
                {
                    breaks.Add(newInd);
                    breaks.Sort();

                    if (verbose == true)
                    {
                        ///
                        //print ("Breakpoint occurs at sample number: ", newInd, ", LL = ",  newVal)
                        // print (len(breaks) - 2, breaks)
                        ///
                    }
                    //Adjust current locations of the breakpoints
                    
                    List<int> lstTemp = new List<int>();
                    lstTemp.Add(newInd);
                    breaks = fAdjustBreaks(data, breaks, lstTemp, lamb, verbose);
                    
                    //Calculate likelihood

                    Double ll = 0.0;
                    ll = fCalculateLikelihood(data, breaks, lamb);
                    breakPoints.Add(breaks);
                    plotPoints.Add(ll);
                }
            }

        }//--------------------------------------------fGgs
        //================================================================================
        //Run cross-validation up to kMax for a set of lambdas
        //Return: train and test set likelihood for every K, lambda

        public void fGGSCrossVal(Double[][] data,int kMax/*=25*/ , Double[] lambList/*=[0.1,1,10]*/ ,
            int[] features/*=[]*/ , bool verbose = false)
        {
            
            if (features != null)
            {
                int nDataLen = data.Length;
                int nFeatureLen = features.Length;
                Double[][] featureData = new Double[nDataLen][];
                for (int ind1 = 0; ind1 < nDataLen; ind1++)
                {
                    featureData[ind1] = new Double[nFeatureLen];
                    for (int ind2 = 0; ind2 < nFeatureLen; ind2++)
                    {
                        featureData[ind1][ind2] = data[ind1][features[ind2]];
                    }
                }
                data = featureData;
            }

            int origSize = data.Length;
            int n = data[0].Length;

            int[] ordering = new int[origSize];
            for(int i=0;i<origSize;i++)
                ordering[i] = i;

            //for each lambda, run the 10 folds in parallel

        }
        //================================================================================
        public List<cMeanCovs> fGGSMeanCov(Double[][] data, List<int> breakpoints, Double lamb, int[] features = null,
            bool verbose=false)
        {
            int m, n;
            List<cMeanCovs> mean_covs = new List<cMeanCovs>();
            int iSize = data[0].Length;

            //Select the desired feature
            if (features != null)
            {
                int nDataLen = data.Length;
                int nFeatureLen = features.Length;
                Double[][] featureData = new Double[nDataLen][];
                for (int ind1 = 0; ind1 < nDataLen; ind1++)
                {
                    featureData[ind1] = new Double[nFeatureLen];
                    for (int ind2 = 0; ind2 < nFeatureLen; ind2++)
                    {
                        featureData[ind1][ind2] = data[ind1][features[ind2]];
                    }
                }
                data = featureData;
            }

            m = data.Length;
            n = data[0].Length;
            int numSegment = breakpoints.Count - 1;

            for (int i = 0; i < numSegment; i++)
            {
                //Get mean and regularized covariance of current segment
                List<Double[]> tempData = new List<double[]>();
                for (int j = breakpoints[i]; j <= breakpoints[i + 1] && j<m; j++)
                    tempData.Add(data[j]);
                m = tempData.Count;
                n = tempData[0].Length;
                ///
                Double[][] dArrayTempData = tempData.ToArray();
                
                Double[] empMean = cStatistic.fMean(dArrayTempData);
                
                Double[][] empCov = cStatistic.fCalcCovariance(dArrayTempData);
                Double[][] mulLambIdenN = cStatistic.fMulScalarMatrix((float)lamb, Matrix.Identity(n).ToJagged());
                mulLambIdenN = cStatistic.fDivMatrixScalar(m, mulLambIdenN);
                Double[][] regularizedCov = cStatistic.fSumTwoMatrix(empCov, mulLambIdenN);
                cMeanCovs oMC = new cMeanCovs();
                oMC.Mean = empMean;
                oMC.RegCovariance = regularizedCov;
                mean_covs.Add(oMC);

            }//end of for (int i = 0; i < numSegment; i++)
            return mean_covs;
        }
        //================================================================================
        public Double fCalculateLikelihood(Double[][] data,List<int> breaks,Double lamb)
        {
            Double ll = 0.0;

            for (int i = 0; i < breaks.Count - 1; i++)
            {
                int nTempDataLen = breaks[i + 1] < data.Length ? breaks[i + 1] - breaks[i] + 1 : data.Length;
                
                List<Double[]> tempData = new List<double[]>();
                for (int j = breaks[i]; j <= breaks[i + 1] && j < data.Length; j++)
                    tempData.Add(data[j]);
                int m = tempData.Count;
                int n = tempData[0].Length;
                Double[][] empCov = cStatistic.fCalcCovariance(tempData.ToArray());
                
                Double[][] mulLambIdenN = cStatistic.fMulScalarMatrix((float)lamb, Matrix.Identity(n).ToJagged());
                mulLambIdenN = cStatistic.fDivMatrixScalar(m, mulLambIdenN);//float(lamb)*np.identity(n)/m
                Double[][] regularizedCov = cStatistic.fSumTwoMatrix(empCov, mulLambIdenN);//empCov + float(lamb)*np.identity(n)/m)
               
                int sign;
                Double slogDetRegCov = cStatistic.fSLogDeterminant(regularizedCov, out sign);//np.linalg.slogdet(empCov + float(lamb)*np.identity(n)/m)

                Double[][] regCovInverse = cStatistic.fInverse(regularizedCov);//np.linalg.inv(empCov + float(lamb)*np.identity(n)/m)

                Double regCovInvTrace = cStatistic.fTrace(regCovInverse);//np.trace(np.linalg.inv(empCov + float(lamb)*np.identity(n)/m))

                ll = ll - (m * slogDetRegCov - (float)lamb * regCovInvTrace);

            }
            return ll;
        }
        //================================================================================

        public void fAddBreak(Double[][] data,Double lamb,out int minIndex, out Double value)
        {
            //Initialize parameters
            int m = data.Length;
            int n = data[0].Length;
            Double[] origMean = cStatistic.fMean(data);
            Double[][] origCov = cStatistic.fCalcCovariance(data);
            Double[][] mulLambIdenN = Matrix.Identity(n).ToJagged();

            for (int index = 0; index < n; index++)
            {
                mulLambIdenN[index][index] = (float)(mulLambIdenN[index][index] * lamb) / m;//float(lamb)*np.identity(n)/m
            }
            Double[][] regularizedCov = cStatistic.fSumTwoMatrix(origCov, mulLambIdenN);//empCov + float(lamb)*np.identity(n)/m)

            int sign;
            Double slogDetRegCov = cStatistic.fSLogDeterminant(regularizedCov, out sign);//np.linalg.slogdet(empCov + float(lamb)*np.identity(n)/m)

            Double[][] regCovInverse = cStatistic.fInverse(regularizedCov);//np.linalg.inv(empCov + float(lamb)*np.identity(n)/m)

            Double regCovInvTrace = cStatistic.fTrace(regCovInverse);//np.trace(np.linalg.inv(empCov + float(lamb)*np.identity(n)/m))
            Double origLL = m * slogDetRegCov - (float)lamb * regCovInvTrace;
            Double[][] origMeanOuter = cStatistic.fOuterProduct(origMean, origMean);
            Double[][] totSum = cStatistic.fSumTwoMatrix(origCov, origMeanOuter);//origCov+np.outer(origMean,origMean)
            totSum = cStatistic.fMulScalarMatrix(m, totSum);//totSum = m*(origCov+np.outer(origMean,origMean))

            Double[] muLeft = cStatistic.fDivArrayScalar(n, data[0]);//muLeft = data[0,:]/n

            Double[] mulMOrigMean = cStatistic.fMulScalarArray(m, origMean);
            Double[] subMulMOmeanData0 = cStatistic.fSubTwoArray(mulMOrigMean, data[0]);//m * origMean - data[0,:]
            Double[] muRight = cStatistic.fDivArrayScalar(m - 1, subMulMOmeanData0);//muRight = (m * origMean - data[0,:])/(m-1)

            Double[][] runSum = cStatistic.fOuterProduct(data[0], data[0]);//runSum = np.outer(data[0,:],data[0,:])
            //
            //Loop through all samples, find point where breaking the segment would have the largest LL increase
            //
            Double minLL = origLL;
            int minInd = 0;
            for (int i = 2; i < m - 1; i++)
            {
                //Update parameters
                runSum = cStatistic.fSumTwoMatrix(runSum, cStatistic.fOuterProduct(data[i - 1], data[i - 1]));
                //runSum = runSum + np.outer(data[i-1,:],data[i-1,:])
                muLeft = cStatistic.fDivArrayScalar(i, cStatistic.fSumTwoArray(cStatistic.fMulScalarArray(i - 1, muLeft), data[i - 1]));
                //muLeft = ((i-1)*muLeft + data[i-1,:])/(i)
                muRight = cStatistic.fDivArrayScalar(m - i, cStatistic.fSubTwoArray(cStatistic.fMulScalarArray(m - i + 1, muRight), data[i - 1]));
                //muRight = ((m-i+1) * muRight - data[i-1,:])/(m-i)
                Double[][] sigLeft = cStatistic.fSubTwoMatrix(cStatistic.fDivMatrixScalar(i, runSum), cStatistic.fOuterProduct(muLeft, muLeft));
                //sigLeft = runSum/(i) - np.outer(muLeft, muLeft)
                Double[][] sigRight = cStatistic.fSubTwoMatrix(cStatistic.fDivMatrixScalar(m - i, cStatistic.fSubTwoMatrix(totSum, runSum)),
                    cStatistic.fOuterProduct(muRight, muRight));
                //sigRight = (totSum - runSum)/(m-i) - np.outer(muRight,muRight)
                //
                //Compute Cholesky, LogDet, and Trace
                //


                Double[][] mulDivIdenForChlskyLleft = cStatistic.fDivMatrixScalar(i, cStatistic.fMulScalarMatrix((float)lamb,
                    Matrix.Identity(n).ToJagged()));//float(lamb)*np.identity(n)/i
                Double[][] Lleft = cStatistic.fCholesky(cStatistic.fSumTwoMatrix(sigLeft, mulDivIdenForChlskyLleft));
                //Lleft = np.linalg.cholesky(sigLeft + float(lamb)*np.identity(n)/i)
                Double[][] mulDivIdenForChlskyLright = cStatistic.fDivMatrixScalar(m - i, cStatistic.fMulScalarMatrix((float)lamb,
                    Matrix.Identity(n).ToJagged()));
                Double[][] Lright = cStatistic.fCholesky(cStatistic.fSumTwoMatrix(sigRight, mulDivIdenForChlskyLright));
                //Lright = np.linalg.cholesky(sigRight + float(lamb)*np.identity(n)/(m-i))

                Double llLeft = 2 * cStatistic.fSumArray(cStatistic.fMapLog(cStatistic.fDiag(Lleft)));//llLeft = 2*sum(map(math.log, np.diag(Lleft)))
                Double llRight = 2 * cStatistic.fSumArray(cStatistic.fMapLog(cStatistic.fDiag(Lright)));//llRight = 2*sum(map(math.log, np.diag(Lright)))

                Double trLeft = 0, trRight = 0;
                if (lamb > 0)
                {
                    trLeft = Math.Pow(cStatistic.fEuclideanNorm(cStatistic.fInverse(Lleft)), 2);
                    trRight = Math.Pow(cStatistic.fEuclideanNorm(cStatistic.fInverse(Lright)), 2);
                }
                Double LL = i * llLeft - (float)lamb * trLeft + (m - i) * llRight - (float)lamb * trRight;
                if (LL < minLL)
                {
                    minLL = LL;
                    minInd = i;
                }
            }
            //Return break, increase in LL
            minIndex = minInd;
            value = minLL - origLL;
        }
        //================================================================================

        public List<int> fAdjustBreaks(Double[][] data, List<int> breakpoints, List<int> newInd, Double lamb = 0.0, bool verbose = false, int maxShuffles = 250)
        {
            List<int> bp = breakpoints;

            if (bp.Count == 3)
                return bp;
            Dictionary<int, int> lastPass = new Dictionary<int, int>();
            Dictionary<int, int> thisPass = new Dictionary<int, int>();
            foreach (int b in bp)//for b in bp
                thisPass[b] = 0;

            foreach (int i in newInd)//for i in newInd
                thisPass[i] = 1;

            for (int z = 0; z < maxShuffles; z++)
            {
                lastPass = thisPass;
                thisPass = new Dictionary<int, int>();
                foreach (int b in bp)
                    thisPass[b] = 0;
                bool switchAny = false;
                List<int> ordering=new List<int>();
                for (int i = 1; i < bp.Count-1; i++)
                    ordering.Add(i);
                ///
                //random.shuffle(ordering)
                ordering.Shuffle();
                ///
                foreach(int i in ordering)
                {
                    //Check if we need to adjust it
                    if((lastPass.ContainsKey(bp[i-1])&&lastPass[bp[i-1]]==1) || lastPass[bp[i+1]]==1
                        || (thisPass.ContainsKey(bp[i - 1])&&thisPass[bp[i-1]]==1) || thisPass[bp[i+1]]==1)
                    {
                        int tempDataLen = bp[i + 1] - bp[i - 1] + 1;

                        List<Double[]> tempData = new List<double[]>();
                        bool b = false;
                        for (int j = bp[i-1]; j <= bp[i + 1] && j < data.Length; j++)
                            tempData.Add(data[j]);
                        int ind; Double val;
                        fAddBreak(tempData.ToArray(), lamb,out ind,out val);
                        if (bp[i] != ind + bp[i - 1] && val != 0)
                        {
                            lastPass[ind + bp[i - 1]] = lastPass[bp[i]];
                            lastPass.Remove(bp[i]);
                            //del lastPass[bp[i]]
                            thisPass.Remove(bp[i]);
                            //del thisPass[bp[i]]

                            thisPass[ind + bp[i - 1]] = 1;
                            bp[i] = ind + bp[i - 1];
                            switchAny = true;
                        }// end of if(bp[i]!=ind+bp[i-1] && val!=0)

                    }//end of if(lastPass[bp[i-1]]==1 || lastPass[bp[i+1]]==1 || thisPass[bp[i-1]]==1 || thisPass[bp[i+1]]==1)
                }//end of foreach(int i in ordering)
                if (switchAny == false)
                    return bp;

            }//end of for (int z = 0; z < maxShuffles; z++)
            return bp;
        }
        //================================================================================
        public void fMulti_run_wrapper(int fold,Double[][] data,int kmax,Double lamb,bool verbose,int origSize,int n,List<int> ordering)
        {
            //******
            
        }
        //================================================================================

        //for identity
        //Matrix.Identity(n);//http://accord-framework.net/docs/html/T_Accord_Math_Matrix.htm
        
       
    }
}
