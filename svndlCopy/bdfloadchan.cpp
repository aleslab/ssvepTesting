//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// bdfloadchan
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "mex.h"
#include "matrix.h"
#include "bdf.h"
//#define	COMPILE_IN_MATLAB
//#include "matlabaux.h"
//#include "matlabaux.cpp"

#ifdef _WIN32
#include <windows.h>

BOOL APIENTRY DllMain( HANDLE hModule, 
                       DWORD  ul_reason_for_call, 
                       LPVOID lpReserved
					 )
{
    return TRUE;
}
#endif

////////////////////////////////////////////////////////////////////////////////
// This is the core of the work
////////////////////////////////////////////////////////////////////////////////
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
	if (nrhs == 1 || nrhs > 3)
	{
		mexEvalString("help bdfloadchan;");
	}
	else
	{
		BdfHdr header;
		
		
		
		//char filename[512];
		char *filename;
		int   buflen;
		
		int chan = maGetINT32Element(prhs[1], 0) - 1;
		int recordsToRead;

		buflen = mxGetN(prhs[0])*sizeof(mxChar)+1;
		filename = (char *)mxMalloc(buflen);
		
		mxGetString(prhs[0], filename, buflen);
		
		
		FILE* fp = fopen(filename, "rb");
		
		if ( fp == NULL )
		{
			printf("ERROR: Unable to load file: %s\n", filename);
		}
		else if (fp && header.Load(fp))
		{
			if (nrhs > 2)
				recordsToRead = maGetINT32Element(prhs[2],0);
			else
				recordsToRead = header.GetNumDataRecord();
		
			int sr;
			int numTotalRecs = header.GetNumDataRecord();
			int nc = header.GetNumChan();
			int i;
			int sampsPerChanPerRecord = header.GetSamplesPerRecord(chan);

			int totalChanSize = sampsPerChanPerRecord * recordsToRead;
			
			
			
		
			
			int adcMax = header.GetAdcMax(chan);
			
			int adcMin = header.GetAdcMin(chan);
			
			
			
			int maxValue = header.GetMaxValue(chan);
			
			
			
			int minValue = header.GetMinValue(chan);
			
			
			
			double conversionFactor = (double)(maxValue-minValue)/(double)(adcMax-adcMin);
			
			double offset = (minValue - (conversionFactor * adcMin));
			
			
			
			char *chanType;
			//chanType = new char[80];
			char chanType2;
			
			chanType = (char *)mxMalloc(81);
			chanType[80] = 0;
			
			header.GetTransducerType(chan,chanType);
			
//			chanType2[0] = 'a';
//			chanType2[1] = 0;
  			
					
			
 			//printf("channel = %d tcs = %d spcpr = %d ntr = %d rtr = %d calib = %f offset =%f\n",
			//		chan, totalChanSize, sampsPerChanPerRecord, 
				//	numTotalRecs, recordsToRead,conversionFactor, offset);
					
					
					
			plhs[0] = mxCreateNumericMatrix(1, totalChanSize, mxDOUBLE_CLASS, mxREAL);
			double* pData = (double*)mxGetData(plhs[0]);
			i = 0;
			// load individual record as they come... into a linked list
			bool done = false;
			char* curRecordData = new char[sampsPerChanPerRecord*3];
			fseek(fp, (nc+1)*256, SEEK_SET);
			for (int k = 0; k < recordsToRead; k++)
			{
				for (int q = 0; q < chan; q++)
					sr = fseek(fp, header.GetSamplesPerRecord(q)*3, SEEK_CUR);
				if (fread(curRecordData, sampsPerChanPerRecord*3, 1, fp))
				{
					unsigned char* temp = (unsigned char*)curRecordData;
					for (int j = 0; j < sampsPerChanPerRecord; j++)
					{
						//unsigned long v = temp[0] | (temp[1] << 8) | (temp[2] << 16);
						int v = temp[0] | (temp[1] << 8) | (temp[2] << 16);
						// convert this from using 24 bits to a full int
						// need to do sign extension
						int r; // resulting sign extended number goes here
						int value;
						value = strncmp(chanType,"Active Electrode",16);
						if (value==0)
						 {
							struct {signed int v:24;} s;
							r = s.v = v;
							pData[i++] = (double)r*conversionFactor + offset;
						}
						else 
						{
						pData[i++] = (double)v;//*conversionFactor + offset;
						}
						
						
						temp += 3;
					}
					for (int q = chan+1; q < nc; q++)
						sr = fseek(fp, header.GetSamplesPerRecord(q)*3, SEEK_CUR);
				}
//				printf("loaded %d records\n", k);
			}
			delete [] curRecordData;
		//	delete [] chanType;
			// all done reading!
		}
		if (fp)
			fclose(fp);

	}
}
