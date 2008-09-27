/*=========================================================================
Program:   FLT File Class
Module:    $RCSfile: FltFile.cpp,v $
Language:  C++
Author:    $Author: bjian $
Date:      $Date: 2007/10/14 17:27:33 $
Version:   $Revision: 1.5 $
=========================================================================*/

#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace std;

#include "myendian.h"
#include "FltFile.h"

//constructor
FltFile::FltFile()
{
    // Try member initialization list
    m_LocalEndian = isThisMachineLittleEndian()?'l':'b';
    m_cReadMode = 0;
    m_cWriteMode = 'b';
    m_bSwap = false;
    m_bOpen = false;
    m_DataType = 4;
    m_nSampTotal = 1;
}

//constructor
FltFile::FltFile(const char* fileName, const char mode)
:m_cReadMode(mode), m_cWriteMode('b'), m_bSwap(false), m_bOpen(false)
{
    // Initialize pointers.
    m_LocalEndian = isThisMachineLittleEndian()?'l':'b';
    m_nSampTotal = 1;
    Read(fileName);
}

//read from file
int FltFile::Read(const char* inputFileName)
{
    ifstream myFile(inputFileName, ios::in | ios::binary);
    if (myFile.is_open())
    {
        if (ReadHeader(myFile))
        {
            // Continue to read data according to different data type.
            switch (m_DataType)
            {
                case DTTYPE_FLOAT:
                    cout << " Datatype is float.\n";
                    return ReadFloat(myFile);
                    break;
                    // case DTTYPE_COMPLEX:
                    //  cout << " Datatype is complex.\n";
                    //  ReadComplex(inputFile, mode);
                    //  break;
                default:
                    printf(" ERROR: Won't read data type %ld.\n", m_DataType);
                    return -1;
                    break;
            }
            return 1;
        }
        else
            return -1;
    }
    else
    {
        cerr << "Unable to open file " << inputFileName << endl;
        return -1;
    }
}

// readheader
int FltFile::ReadHeader(ifstream& myFile)
{
    unsigned dim;
    bool swapflg = 0;
    myFile.read((char *)&dim, 4); // 4==sizeof(unsigned)
    m_bSwap =(dim>IMPOSSIBLY_LARGE_DIM);
    if (m_bSwap)
    {
        swapbytes((unsigned char*)&dim, 4);
    }
    m_nDimTotal = dim;
    if (m_nDimTotal != 4 && m_nDimTotal != 3)
    {
        cerr << " ERROR: Can't read non-4D or non-3D data set.\n";
        return -1;
    }
    unsigned size_, nElemsTotal =1;
    m_nVolume.push_back(nElemsTotal);
    for (unsigned i = 0; i < m_nDimTotal; ++i)
    {
        myFile.read((char*)&size_, 4);
        if (m_bSwap)
            swapbytes((unsigned char*)&size_, 4);
        nElemsTotal *= size_;
        m_nSize.push_back(size_);
        m_nVolume.push_back(nElemsTotal);
    }
    if (m_nDimTotal == 4)
    {
        //For 3D data, sample total is set to 1 in the constructor
        m_nSampTotal = m_nSize[3];
    }
    myFile.read((char*) &m_DataType, 4);
    myFile.read((char*) &m_nElemsTotal, 4);
    if (m_bSwap)
    {
        swapbytes((unsigned char*)&m_nElemsTotal, 4);
        swapbytes((unsigned char*)&m_DataType, 4);
    }
    return 1;
}

//readfloat
int FltFile::ReadFloat(ifstream& myFile)
{
    if (m_DataType != DTTYPE_FLOAT)
    {
        cerr << " ERROR: Data type is not 32bit float.\n";
        return -1;
    }
    for (unsigned i = 0; i < m_nSampTotal; ++i)
    {
        float* m_pVolume = new float[m_nVolume[3]];
        myFile.read((char*)m_pVolume, 4*m_nVolume[3]);
        if (m_cReadMode == 0 && m_bSwap)
        {
            //default option is to assume the byteorder in data section is the same as in the header
            for (unsigned j = 0; j < m_nVolume[3]; ++j)
            {
                swapbytes((unsigned char*)(m_pVolume + j), 4);
            }
        }
        if (m_cReadMode>0 && m_cReadMode != m_LocalEndian)
        {
            //if the read mode is explicitly specified, check if it is consistent with the local endian
            for (unsigned j = 0; j < m_nVolume[3]; ++j)
            {
                swapbytes((unsigned char*)(m_pVolume + j), 4);
            }
        }
        m_Data.push_back(m_pVolume);
    }
    m_bOpen = true;
    return 1;
}

//readsignal // position starts from 0
int FltFile::ReadSignal(std::vector<float>& signal, unsigned position)
{
    signal.clear();
    if (position >= 0 && position < m_nVolume[3])
    {
        for (unsigned i = 0; i < m_nSampTotal; ++i)
        {
            signal.push_back(*((m_Data[i]) + position));
        }
    }
    else
    {
        cerr << "position out of range\n";
        return -1;
    }
    return 1;
}

//readsignal;  //coordinate starts from [0,0,0]
int FltFile::ReadSignal(std::vector<float>& signal, unsigned* coordinate)
{
    unsigned position = 0;

    for (unsigned i = 0; i < 3; ++i) //assume the dimension of volume is always 3
    {
        if (coordinate[i] >= 0 && coordinate[i] < m_nSize[i])
            position += coordinate[i]*m_nVolume[i];
        else
        {
            cerr << "coordinate out of range\n";
            return -1;
        }
    }
    return ReadSignal(signal, position);
}

// Destructor.
FltFile::~FltFile()
{
    emptyHeader();
}

// Free allocated memory.
void FltFile::emptyHeader()
{
    /*
    for (unsigned i = 0; i < m_Data.size(); ++i)
    {
        delete[] m_Data[i];
    }*/
    for (std::vector<float*>::iterator pObj=m_Data.begin(); pObj!=m_Data.end(); ++pObj)
    {
        delete[] *pObj;
    }
    m_Data.clear();
    m_nSize.clear();
    m_nVolume.clear();
}


// Write flt data to a file.
int FltFile::WriteFile(const char* outputFileName, std::vector<float*>& data_, std::vector < unsigned>& size_, const char mode, bool b2DSlice)
{
    ofstream outFile(outputFileName, ios::out | ios::binary);
    if (outFile.is_open())
    {
        emptyHeader();
        m_nSampTotal = 1;
        //this->m_Data = data_;
        this->m_nSize = size_;
        m_nDimTotal = m_nSize.size();
        // check integrity and other numeric values!!
        unsigned nElemsTotal = 1;
        m_nVolume.push_back(nElemsTotal);

        for (unsigned i = 0; i < m_nDimTotal; ++i)
        {
            nElemsTotal *= m_nSize[i];
            m_nVolume.push_back(nElemsTotal);
        }

        std::cout << "m_nVolume: [";
        for (unsigned k=0; k<m_nVolume.size(); ++k)
        {
            std::cout << m_nVolume[k] << " ";
        }
        std::cout << "]"<<std::endl;

        m_nElemsTotal = nElemsTotal;
        if (m_nDimTotal == 4)
        {
             m_nSampTotal = m_nSize[3]; //default is 1 for 3D.
        }
        if (m_nDimTotal ==3 && b2DSlice )
        {
            m_nSampTotal = m_nSize[2]; // now assume 2D slice and 3rd dim is nSampTotal
        }
//      if (m_nDimTotal ==3 && !b2DSlice )
        {
        //  m_nSampTotal = 1; // now assume 3D slice and 3rd dim is nSampTotal
        }

        std::cout << "m_nElemsTotal "<< m_nElemsTotal<<std::endl;
        std::cout << "m_nSampTotal "<< m_nSampTotal<<std::endl;

        m_cWriteMode = mode;
        m_bOpen = true;
        // ! debug only: PrintHeader();

        // Write the header.
        if (WriteHeader(outFile))
        {
            // Continue to write data according to different data type.
            switch (m_DataType)
            {
                case DTTYPE_FLOAT:
                    WriteFloat(outFile,data_);
                    break;
                /* Implement later
                case DTTYPE_COMPLEX:
                WriteComplex(outputFile);
                break;
                */
                default:
                    cerr<<" ERROR: Won't read data type " << m_DataType << endl;
                    return -1;
                    break;
            }
        }
        else
        {
            cerr<<" Errors occured when trying to write the header\n ";
            return -1;
        }
        return 1;
    }
    else
    {
        cerr << "Unable to open file " << outputFileName << endl;
        return -1;
    }
}

//write header
int FltFile::WriteHeader(ofstream& outFile)
{
    // Write the total number of dimensions from file
    // Only handle 4D data now
    if (m_nSize.size() != 4 && m_nSize.size() != 3)
    {
        cerr << " ERROR: Can't write non-4D or non-3D data set.\n";
        return (-1);
    }
    unsigned temp = m_nDimTotal;
    if (m_cWriteMode != m_LocalEndian)
        swapbytes((unsigned char*)&temp, 4);
    outFile.write((char*)&temp, 4);

    for (unsigned i = 0; i < m_nDimTotal; ++i)
    {
        temp = m_nSize[i];
        if (m_cWriteMode != m_LocalEndian)
            swapbytes((unsigned char*)&temp, 4);
        outFile.write((char*)&temp, 4);
    }
    temp = m_DataType;
    if (m_cWriteMode != m_LocalEndian)
    {
        swapbytes((unsigned char*)&temp, 4);
    }
    outFile.write((char*)&temp, 4);
    temp = m_nElemsTotal;
    if (m_cWriteMode != m_LocalEndian)
    {
        swapbytes((unsigned char*)&temp, 4);
    }
    outFile.write((char*)&temp, 4);
    return 1;
}

//write float
void FltFile::WriteFloat(ofstream& outFile)
{
    WriteFloat(outFile, m_Data);
}

void FltFile::WriteFloat(ofstream& outFile, vector<float*>& data)
{
    for (unsigned i = 0; i < m_nSampTotal; ++i)
    {
        float* m_pVolume = data[i];
        for (unsigned j = 0; j < m_nVolume[m_nDimTotal]/m_nSampTotal; ++j)
        {
            float temp = m_pVolume[j];
            if (m_cWriteMode != m_LocalEndian)
                 swapbytes((unsigned char*)&temp, 4);
            outFile.write((char*)&temp, 4);
        }
    }
}

//print header
void FltFile::PrintHeader()
{
    if (m_bOpen)
    {
        cout << "=====================================================\n";
        cout << " Total number of dimensions : " << m_nDimTotal << endl;
        unsigned nElemsTotal = 1;
        for (unsigned i = 0; i < m_nSize.size(); ++i)
        {
            cout << " Size in dim[" << i << "]: " <<  m_nSize[i] << endl;
            nElemsTotal *= m_nSize[i];
        }
        cout << " Calculated total number of data elements is " << nElemsTotal << endl;
        cout << " Data type is " << m_DataType << endl;
        cout << " Total number of data elements is " << m_nElemsTotal << endl;
        cout << " # of samples is " << m_nSampTotal << endl;
        for (unsigned j = 0; j < m_nVolume.size(); ++j)
        {
            cout << " Number in m_nVolume[" << j << "]: " <<  m_nVolume[j] << endl;
        }
        cout << "=====================================================\n";

    }
    else
    {
        cerr << "Empty Object!\n";
    }
}


void FltFile::CheckData()
{
    if (is_open())
    {
        std::vector<float> signal;
        unsigned coord[3] = {1,0,0};
        ReadSignal(signal, coord);
        std::cout<<"The signal at second voxel [1,0,0]:" <<std::endl;
        for (unsigned int k=0;k<signal.size();++k)
            std::cout<<signal[k]<<" ";
        std::cout<<std::endl;
        unsigned position = volume()/2;
        ReadSignal(signal, position);
        std::cout<<"The signal at position ("<< position<<"):\n";
        for (unsigned int k=0;k<signal.size();++k)
            std::cout<<signal[k]<<" ";
        std::cout<<std::endl;
        position = volume()-1;
        ReadSignal(signal, position);
        std::cout<<"The signal at last voxel ("<< position<<"):\n";
        for (unsigned int k=0;k<signal.size();++k)
            std::cout<<signal[k]<<" ";
        std::cout<<std::endl;
    }
}
