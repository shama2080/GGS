using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GGS
{
    public class cMeanCovs
    {
        Double[] _dMean;
        Double[][] _dRegCovariance;

        //--------------------------------------------------------------------------------
        public Double[] Mean 
        {
            set
            {
                this._dMean = value;
            }
            get
            {
                return this._dMean;
            }
        }

        //--------------------------------------------------------------------------------
        public Double[][] RegCovariance
        {
            set
            {
                this._dRegCovariance = value;
            }
            get
            {
                return this._dRegCovariance;
            }
        }
        //--------------------------------------------------------------------------------
    }
}
