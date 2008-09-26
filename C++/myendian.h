/*=========================================================================
Program:   FLT File Class
Module:    $RCSfile: myendian.h,v $
Language:  C
Author:    $Author: bjian $
Date:      $Date: 2006/09/21 19:53:59 $
Version:   $Revision: 1.1 $
=========================================================================*/

#ifndef _MY_ENDIAN_H_
#define _MY_ENDIAN_H_

# ifdef __cplusplus
  extern "C"{
# endif


int isThisMachineLittleEndian();
void swapbytes(unsigned char*, unsigned size);

# ifdef __cplusplus
  }
# endif

#endif
/* End of myendian.h */
