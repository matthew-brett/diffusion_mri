/*=========================================================================
Program:   FLT File Class
Module:    $RCSfile: FltFile.h,v $
Language:  C++
Author:    $Author: bjian $
Date:      $Date: 2007/08/10 03:44:19 $
Version:   $Revision: 1.3 $
=========================================================================*/

// Based on code written by Tim McGraw and Zhizhou Wang

#ifndef _FLT_FILE_H__
#define _FLT_FILE_H__

#include <iostream>
#include <fstream>
#include <vector>

#define _LITTLE_ENDIAN
#define LITTLE_ENDIAN_MODE		0
#define BIG_ENDIAN_MODE			1

const unsigned IMPOSSIBLY_LARGE_DIM = 100;
typedef unsigned char BYTE;
typedef unsigned int  WORD32;

enum {DTTYPE_UNDEFINED, DTTYPE_BYTE, DTTYPE_INTEGER, DTTYPE_LWORD, DTTYPE_FLOAT, DTTYPE_DOUBLE, DTTYPE_COMPLEX};
/*
 * 7       String
 * 8       Structure
 * 9       Double-precision complex
 * 10      Pointer
 * 11      Object reference
 * 12      Unsigned Integer
 * 13      Unsigned Longword Integer
*/


class FltFile
{
	protected:
		// for file header
		unsigned m_nDimTotal;
  		unsigned m_DataType;
  		unsigned m_nElemsTotal;
		unsigned m_nSampTotal;
		std::vector<float *> m_Data;
		std::vector<unsigned> m_nSize;
		std::vector<unsigned> m_nVolume;

	private:
		int Read(const char*);
		int ReadHeader(std::ifstream&);
		int ReadFloat(std::ifstream&);
		int WriteHeader(std::ofstream&);
		void WriteFloat(std::ofstream&);
		void emptyHeader();
		void WriteFloat(std::ofstream&, std::vector<float*>&);

		mutable char m_LocalEndian, m_cReadMode, m_cWriteMode;
		mutable bool m_bSwap, m_bOpen;
	public:
  		FltFile();
		FltFile(const char*, const char = 0);
		int WriteFile(const char*, std::vector<float*>&, std::vector<unsigned>&, const char = 'b',bool b2DSlice=0);
  		~FltFile();
		
  		// Read, peek and write data from flt file;
  		/*int Read(char* fileName);
  		int Read(char* fileName, int mode);
  		int Peek(char* fileName);
  		int Write(char* fileName);
  		int Write(char* fileName, int mode);

		WORD32 ReadFromMemory(float* data, int* nsize);
		WORD32 ReadFromMemory(float* data[7], int* nsize);
	
  		void GetSize(int size[4]);
  		float* GetDataMag(int nSampId);
		float* GetDataMag();
  		float* GetDataPhase(int nSampId);
  		float* GetDataReal(int nSampId);
  		float* GetDataImg(int nSampId);

  		int GetSampTotal(void) 
		{
			return m_nSampTotal;
		}
  		inline int Index (int i, int j, int k, int l) {
			return ((l*m_nDepth+k)*m_nHeight+j)*m_nWidth + i;
		};
  		inline int Index3D (int i, int j, int k) {
			return (k*m_nHeight+j)*m_nWidth + i;
		};
		WORD32 GetDataType();
		int m_SourceHeadEndianMode;
		int m_SourceDataEndianMode;
		int m_ThisMachineEndianMode;
		int m_TargetHeadEndianMode;
		int m_TargetDataEndianMode;
		void SetSourceDataEndianMode(int mode);
		void SetTargetDataEndianMode(int mode);
        */

		inline bool is_open(){
			return m_bOpen;
		}
  		inline int samples(){
			return m_nSampTotal;
		}
		inline std::vector<float*> & getData(){
			return m_Data;
		}
		inline std::vector<unsigned> & getSize(){
			return m_nSize;
		}
		inline int getSamples(){
			return m_nSampTotal;
		}
		inline int getVolume(){
			if (m_bOpen)
				return m_nVolume[3];
			else
				return -1;	
		}
  		inline int volume(){
			return getVolume();
		}
		/*
		inline void setData(std::vector<float*>& data_){
			//delete old;
			this->m_Data = data_;
		}
		inline void setSize(std::vector<unsigned>& size_){
			//delete old;
			this->m_nSize = size_;
		}*/

		void CheckData();
		int ReadSignal(std::vector<float>&, unsigned);
		int ReadSignal(std::vector<float>&, unsigned*);
		void PrintHeader();
};

#endif
