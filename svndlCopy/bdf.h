#ifndef BDF_H
#define BDF_H

#include <stdio.h>
#include <string.h>
#include <memory.h>
#include <stdlib.h>
#ifdef MEX_COMPILE
#include "mex.h"
#include "matrix.h"
#endif

#ifdef _MSC_VER
#pragma pack(push, 1)
#endif

struct BdfInfo
{
	char magic[8];
	char subject[80];
	char recording[80];
	char startdate[8];
	char starttime[8];
	char headerSize[8];
	char formatVersion[44];
	char numDataRecord[8];
	char dataDuration[8];
	char numChan[4];
};

struct BdfHdr
{
	BdfInfo info;

	// information about the channels
	char* label;
	char* transducerType;
	char* dimension;
	char* minValue;
	char* maxValue;
	char* adcMin;
	char* adcMax;
	char* prefiltering;
	char* samplesPerRecord;

	BdfHdr()
	{
		label = 
		transducerType = 
		dimension = 
		minValue = 
		maxValue = 
		adcMin = 
		adcMax = 
		prefiltering = 
		samplesPerRecord = 0;
		memset(&info, 0, sizeof(info));
	}

	~BdfHdr()
	{
		if (label) delete [] label;
		if (transducerType) delete [] transducerType;
		if (dimension) delete [] dimension;
		if (minValue) delete [] minValue;
		if (maxValue) delete [] maxValue;
		if (adcMin) delete [] adcMin;
		if (adcMax) delete [] adcMax;
		if (prefiltering) delete [] prefiltering;
		if (samplesPerRecord) delete [] samplesPerRecord;
	}

	int GetSamplesPerRecord(int chan)
	{
		char temp[9]; temp[8] = 0;
		memcpy(temp,&samplesPerRecord[8*chan],sizeof(temp)-1);
		return atoi(temp);
	}
	
	int GetNumChan()
	{
		char temp[5]; temp[4] = 0;
		memcpy(temp, info.numChan, sizeof(temp)-1);
		return atoi(temp);
	}

	int GetAdcMin(int chan)
	{
		char temp[9]; temp[8] = 0;
		memcpy(temp, &adcMin[8*chan], sizeof(temp)-1);
		return atoi(temp);
	}

	int GetAdcMax(int chan)
	{
		char temp[9]; temp[8] = 0;
		memcpy(temp, &adcMax[8*chan], sizeof(temp)-1);
		return atoi(temp);
	}
	
	int GetMaxValue(int chan)
	{
		char temp[9]; temp[8] = 0;
		memcpy(temp, &maxValue[8*chan], sizeof(temp)-1);
		return atoi(temp);
	}

	int GetMinValue(int chan)
	{
		char temp[9]; temp[8] = 0;
		memcpy(temp, &minValue[8*chan], sizeof(temp)-1);
		return atoi(temp);
	}
	
	
	
	void GetTransducerType(int chan, char* chanType)
	{
		
		//char temp[81];temp[80] = 0;
		
		
		strncpy(chanType,&transducerType[80*chan],80);
		//memcpy(temp,&transducerType[0],sizeof(temp)-1);
		
		
		return;
	}



	int GetNumDataRecord()
	{
		char temp[9]; temp[8] = 0;
		memcpy(temp,info.numDataRecord,sizeof(temp)-1);
		return atoi(temp);
	}
	
	
	void Alloc(int nChan)
	{
		// allocate data appropriately
		label = new char[nChan * 16];
		transducerType = new char[nChan * 80];
		dimension = new char[nChan * 8];
		minValue = new char[nChan * 8];
		maxValue = new char[nChan * 8];
		adcMin = new char[nChan * 8];
		adcMax = new char[nChan * 8];
		prefiltering = new char[nChan * 80];
		samplesPerRecord = new char[nChan * 8];
	}

	//#ifdef MEX_COMPILE
	mxArray* ToMatlabArray()
	{
		mxArray* pResult = 0;

	}
	//#endif


	bool Load(FILE* fp)
	{
		bool bResult = false;
		if (fp)
		{
			if (fread(&info, sizeof(info), 1, fp))
			{
				int nChan = atoi(info.numChan);

				Alloc(nChan);
				if (fread(label, nChan * 16, 1, fp) &&
				fread(transducerType, nChan * 80, 1, fp) &&
				fread(dimension, nChan * 8, 1, fp) &&
				fread(minValue, nChan * 8, 1, fp) &&
				fread(maxValue, nChan * 8, 1, fp) &&
				fread(adcMin, nChan * 8, 1, fp) &&
				fread(adcMax, nChan * 8, 1, fp) &&
				fread(prefiltering, nChan * 80, 1, fp) &&
				fread(samplesPerRecord, nChan * 8, 1, fp))
				{
					char junk[32];
					bResult = true;
					for (int i = 0; i < nChan; i++)
					{
						if (!fread(junk, 32, 1, fp))
						{
							bResult = false;
							break;
						}
					}
				}
			}
		}
		return bResult;
	}
};


#ifdef _MSC_VER
#pragma pack(pop)
#endif

#endif
