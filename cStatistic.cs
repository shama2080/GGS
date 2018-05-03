using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Accord.Statistics;
using Accord.Math;
using Accord.Math.Decompositions;

namespace GGS
{
    public static class cStatistic
    {
        //================================================================================
        //---------------------------------Covariance-------------------------------------
        public static Double fCalcCovariance(Double[] dArray)
        {
            return dArray.Covariance(dArray,false);//Accord.Statistics (unbiased=false is equal biase=true in python) 
        }
        //================================================================================
        public static Double fCalcCovariance(Double[] dArray1,Double[] dArray2)
        {
            return dArray1.Covariance(dArray2,false);//Accord.Statistics
        }
        //================================================================================
        public static Double[][] fCalcCovariance(Double[][] dMatrix)
        {

            int nMatrixCol = dMatrix[0].Length;
            int nMatrixRow = dMatrix.Length;
            Double[][] dTempMatrix = new Double[nMatrixCol][];
            dTempMatrix = Matrix.Transpose(dMatrix);
            Double[][] dResult = new Double[nMatrixCol][];
            for (int i = 0; i < nMatrixCol; i++)
            {
                dResult[i] = new Double[nMatrixCol];
                for (int j = 0; j < nMatrixCol; j++)
                {
                    dResult[i][j] = (float)dTempMatrix[i].Covariance(dTempMatrix[j], false);
                }
            }
            return dResult;
        }
        //================================================================================
        public static Double[][] fCalcCovariance(Double[][] dMatrix1,Double[][] dMatrix2)
        {
            int nMatrix1Len = dMatrix1.Length;
            int nMatrix2Len = dMatrix2.Length;
            if(nMatrix1Len!=nMatrix2Len)
            {
                throw new Exception("Dimentions of your matrixes is not valid!!!");
            }
            Double[][] dResult = new Double[nMatrix1Len][];
            for (int i = 0; i < nMatrix1Len; i++)
            {
                dResult[i] = new Double[nMatrix1Len];
                for (int j = 0; j < nMatrix1Len; j++)
                {
                    dResult[i][j] = dMatrix1[i].Covariance(dMatrix2[j],false);//Accord.Statistics
                }
            }
            return dResult;
        }
        //================================================================================
        //--------------------------------------Mean with axis=0--------------------------
        public static Double[] fMean(Double[][] dMatrix)//My function. This function calculates mean on cols
        {
            int dRowLen = dMatrix.Length;
            int dColLen = dMatrix[0].Length;
            Double[] dResult = new Double[dColLen];
            Double dSum;
            for (int j = 0; j < dColLen; j++)
            {
                dSum = 0.0;
                for (int i = 0; i < dRowLen; i++)
                {
                    dSum = dSum + dMatrix[i][j];
                }
                dResult[j] = dSum / dRowLen;
            }
            return dResult;
        }
        //================================================================================
        //--------------------------------------Cholesky----------------------------------
        public static double[][] fCholesky(double[][] dMatrix)
        {
            //From https://rosettacode.org/wiki/Cholesky_decomposition
            int n = dMatrix.Length;

            double[][] ret = new double[n][];
            for (int r = 0; r < n; r++)
            {
                ret[r] = new double[n];
                for (int c = 0; c <= r; c++)
                {
                    if (c == r)
                    {
                        double sum = 0;
                        for (int j = 0; j < c; j++)
                        {
                            sum += ret[c][j] * ret[c][j];
                        }
                        ret[c][c] = Math.Sqrt(dMatrix[c][c] - sum);
                    }
                    else
                    {
                        double sum = 0;
                        for (int j = 0; j < c; j++)
                            sum += ret[r][j] * ret[c][j];
                        ret[r][c] = 1.0 / ret[c][c] * (dMatrix[r][c] - sum);
                    }
                }
            }

            return ret;
           
        }
        //================================================================================
        //--------------------------------------Inverse of a Matrix-----------------------
        public static Double[][] fInverse(Double[][] dMatrix)
        {
            Double[][] dResult;
            dResult = Matrix.Inverse(dMatrix);//Accord.Math
            return dResult;
        }
        //================================================================================
        //----------------------------------------Trace of a Matrix-----------------------
        public static Double fTrace(Double[][] dMatrix)//My function
        {
            //Return the sum along diagonal of the array.
            Double dResult = 0.0;
            int nMtxLen = dMatrix.Length;
            for (int i = 0; i < nMtxLen; i++)
            {
                dResult = dResult + dMatrix[i][i];
            }
            return dResult;
        }
        //================================================================================
        //----------------------------------------norm of a Matrix------------------------
        public static Double fEuclideanNorm(Double[][] dMatrix)
        {
            Double dResult = Norm.Euclidean(dMatrix);//Accord.Math
            return dResult;
        }
        //================================================================================
        //-------------------------------------Quick Sort Algorithm-----------------------
        public static void fQuickSort(Double[][] dMatrix)//My function
        {
            int nMatrixRowsLen = dMatrix.Length;
            int nMatrixColsLen = dMatrix[0].Length;
            for (int i = 0; i < nMatrixRowsLen; i++)
            {
                fQuickSort(dMatrix[i], 0, nMatrixColsLen-1);
            }
        }
        //================================================================================
        private static void fQuickSort(Double[] arr, int start, int end)
        {
        //http://csharpexamples.com/c-quick-sort-algorithm-implementation/
            int i;
            if (start < end)
            {
                i = fPartition(arr, start, end);

                fQuickSort(arr, start, i - 1);
                fQuickSort(arr, i + 1, end);
            }
        }
        //================================================================================
        private static int fPartition(Double[] arr, int start, int end)
        {
            Double temp;
            Double p = arr[end];
            int i = start - 1;

            for (int j = start; j <= end - 1; j++)
            {
                if (arr[j] <= p)
                {
                    i++;
                    temp = arr[i];
                    arr[i] = arr[j];
                    arr[j] = temp;
                }
            }

            temp = arr[i + 1];
            arr[i + 1] = arr[end];
            arr[end] = temp;
            return i + 1;
        }
        //================================================================================
        //-------------------------------------Outer Product------------------------------

        public static Double[][] fOuterProduct(Double[] dArray1,Double[] dArray2)
        {
            Double[,] dTempResult = Matrix.Outer(dArray1, dArray2);//Accord.Math
            return dTempResult.ToJagged();
        }
        //================================================================================
        //----------------------------------Ones(dType=bool)------------------------------

        public static bool[][] fOnes_bool(int size)//My function
        {
            bool[][] bResult = new bool[size][];
            for (int i = 0; i < size; i++)
            {
                bResult[i] = new bool[size];
                for (int j = 0; j < size; j++)
                {
                    bResult[i][j] = true;
                }
            }
            return bResult;
        }
        //================================================================================
        //--------------Diag(return a diagonal or construct a diagonal array--------------

        public static Double[] fDiag(Double[][] dMatrix)//My function
        {
            int nRowSize = dMatrix.Length;
            int nColSize = dMatrix[0].Length;
            if (nRowSize != nColSize) throw new Exception("Enter Square matrix!!");
            Double[] dResult = new Double[nRowSize];
            for (int i = 0; i < nRowSize; i++)
            {
                dResult[i] = dMatrix[i][i];
            }
            return dResult;
        }
        //================================================================================
        //------------------------------------Log Determinant-----------------------------
        public static Double fSLogDeterminant(Double[][] dMatrix, out int sign)
        {
            Double[,] ddMatrix = dMatrix.ToMatrix();
            Double dDeterminant = Matrix.Determinant(ddMatrix);
            sign = dDeterminant < 0 ? -1 : +1;
            Double dLogDeterminant = Math.Log(Math.Abs(dDeterminant));
            return dLogDeterminant;
        }

        //================================================================================
        //------------------------------------Random.Shuffle------------------------------


        //----In python we have 'np.random.seed(0)' and 'random.shuffle(ordering)'

        static Random _random = new Random(0);
        static Double[] RandomizeStrings(Double[] arr)
        {
            //https://www.dotnetperls.com/shuffle
            List<KeyValuePair<int, Double>> list =
                new List<KeyValuePair<int, Double>>();
            // Add all strings from array.
            // ... Add new random int each time.
            foreach (Double d in arr)
            {
                list.Add(new KeyValuePair<int, Double>(_random.Next(), d));
            }
            // Sort the list by the random number.
            var sorted = from item in list
                         orderby item.Key
                         select item;
            // Allocate new string array.
            Double[] result = new Double[arr.Length];
            // Copy values to array.
            int index = 0;
            foreach (KeyValuePair<int, Double> pair in sorted)
            {
                result[index] = pair.Value;
                index++;
            }
            // Return copied array.
            return result;
        }

        //================================================================================
        //---------------------------------Sum of two matrix------------------------------
        public static Double[][] fSumTwoMatrix(Double[][] dMatrix1,Double[][] dMatrix2)//My Function
        {
            //This function is our need and was not in python version
            int nRowLen1 = dMatrix1.Length;
            int nRowLen2 = dMatrix2.Length;
            int nColLen1 = dMatrix1[0].Length;
            int nColLen2 = dMatrix2[0].Length;
            if (nColLen1 != nColLen2 || nRowLen1 != nRowLen2)
                throw new Exception("Dimensions of two matrix should be equal");
            Double[][] dSum = new Double[nRowLen1][];
            for (int i = 0; i < nRowLen1; i++)
            {
                dSum[i] = new Double[nColLen1];
                for (int j = 0; j < nColLen1; j++)
                {
                    dSum[i][j] = dMatrix1[i][j] + dMatrix2[i][j];
                }
            }
            return dSum;
        }
        //================================================================================
        //---------------------------------Sum of two matrix------------------------------
        public static Double[] fSumTwoArray(Double[] dArray1, Double[] dArray2)//My Function
        {
            //This function is our need and was not in python version
            int nRowLen1 = dArray1.Length;
            int nRowLen2 = dArray2.Length;
            if (nRowLen1 != nRowLen2)
                throw new Exception("Dimensions of two Array should be equal");
            Double[] dSum = new Double[nRowLen1];
            for (int i = 0; i < nRowLen1; i++)
            {
                dSum[i] = dArray1[i] + dArray2[i];
            }
            return dSum;
        }
        //================================================================================
        //-------------------------Multiply scalar to matrix------------------------------
        public static Double[][] fMulScalarMatrix(Double nScalar,Double[][] dMatrix)
        {
            int nRow = dMatrix.Length;
            int nCol = dMatrix[0].Length;
            Double[][] dResult = new Double[nRow][];
            for (int i = 0; i < nRow; i++)
            {
                dResult[i] = new Double[nCol];
                for (int j = 0; j < nCol; j++)
                {
                    dResult[i][j] = nScalar * dMatrix[i][j];
                }
            }
            return dResult;
        }
        //================================================================================
        //-------------------------Multiply scalar to Array-------------------------------
        public static Double[] fMulScalarArray(Double nScalar, Double[] dArray)
        {
            int nRow = dArray.Length;
            Double[] dResult = new Double[nRow];
            for (int i = 0; i < nRow; i++)
            {
                dResult[i] = nScalar * dArray[i];
                
            }
            return dResult;
        }
        //================================================================================
        //------------------------------Subtract two mutrix--------------------------------
        public static Double[][] fSubTwoMatrix(Double[][] dMatrix1,Double[][] dMatrix2)
        {
            int nRow = dMatrix1.Length;
            int nCol = dMatrix1[0].Length;
            Double[][] dResult = new Double[nRow][];
            for (int i = 0; i < nRow; i++)
            {
                dResult[i] = new Double[nCol];
                for (int j = 0; j < nCol; j++)
                {
                    dResult[i][j] = dMatrix1[i][j] - dMatrix2[i][j];
                }

            }
            return dResult;
        }
        //================================================================================
        //---------------------------Subtract two Array-----------------------------------
        public static Double[] fSubTwoArray(Double[] dArray1, Double[] dArray2)
        {
            int nRow = dArray1.Length;
            Double[] dResult = new Double[nRow];
            for (int i = 0; i < nRow; i++)
            {
                dResult[i] = dArray1[i] - dArray2[i];

            }
            return dResult;
        }
        //================================================================================
        //-------------------------Divide Matrix by scalar--------------------------------
        public static Double[][] fDivMatrixScalar(int nScalar, Double[][] dMatrix)
        {
            int nRow = dMatrix.Length;
            int nCol = dMatrix[0].Length;
            Double[][] dResult = new Double[nRow][];
            for (int i = 0; i < nRow; i++)
            {
                dResult[i] = new Double[nCol];
                for (int j = 0; j < nCol; j++)
                {
                    dResult[i][j] = dMatrix[i][j]/nScalar;
                }

            }
            return dResult;
        }
        //================================================================================
        //-------------------------Divide Array by scalar---------------------------------
        public static Double[] fDivArrayScalar(int nScalar,Double[] dArray)
        {
            int nRow = dArray.Length;
            Double[] dResult = new Double[nRow];
            for (int i = 0; i < nRow; i++)
            {
                dResult[i] = dArray[i]/nScalar;

            }
            return dResult;
        }
        //================================================================================
        //-------------------------Map Log To All Elements--------------------------------
        public static Double[] fMapLog(Double[] dArray)
        {
            int nRow = dArray.Length;
            Double[] dResult = new Double[nRow];
            for(int i=0;i<nRow;i++)
            {
                dResult[i] = Math.Log(dArray[i]);
            }
            return dResult;
        }
        //================================================================================
        //-----------------------------Sum Of Array Elements------------------------------
        public static Double fSumArray(Double[] dArray)
        {
            Double dSum = 0.0;
            foreach (Double d in dArray)
            {
                dSum = dSum + d;
            }
            return dSum;
        }
        //================================================================================
        //----------------------The data should be in column form-------------------------
        public static double[][] fReshape(double[][] dData)
        {
            //double[][] dReshapedData;
            int nRow = dData.Length;
            int nCol = dData[0].Length;
            if (nRow<nCol)
            {
                dData = Matrix.Transpose(dData);
            }
            return dData;
        }
    }
}
