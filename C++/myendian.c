/*=========================================================================
Program:   FLT File Class
Module:    $RCSfile: myendian.c,v $
Language:  C
Author:    $Author: bjian $
Date:      $Date: 2006/09/21 19:53:44 $
Version:   $Revision: 1.1 $
=========================================================================*/

#include "myendian.h"

/*
* isThisMachineLittleEndian
*
* Function attempts to determine if this machine stores
* its values in little endian notation or not.  This is
* used to determine whether the MIPS data must be
* reversed or not.
*
* Function determines endianness by storing a known
* sequence of bits and then observing what value is at
* the "left" end.
*/
int isThisMachineLittleEndian()
{
	int x = 0x87654321;
	char* t =(char*) &x;
	return (t[0] ==(char)0x21);
}



// Use to solve the big-endian and little-endian problem.
void swapbytes(unsigned char*x, unsigned size)
{
	unsigned char c;
	unsigned short s;
	unsigned int l;
	switch (size)
	{
    case 2: /* swap two bytes */
		c = *x;
		*x = *(x + 1);
		*(x + 1) = c;
		break;
    case 4: /* swap two shorts(2 - byte words) */
		s = *(unsigned short *)x;
		*(unsigned short *)x = *((unsigned short *)x + 1);
		*((unsigned short *)x + 1) = s;
		swapbytes((unsigned char *)x, 2);
		swapbytes((unsigned char *)((unsigned short *)x + 1), 2);
		break;
    case 8: /* swap two ints(4 - bytes words) */
		l = *(unsigned int *)x;
		*(unsigned int *)x = *((unsigned int *)x + 1);
		*((unsigned int *)x + 1) = l;
		swapbytes((unsigned char *)x, 4);
		swapbytes((unsigned char *)((unsigned int *)x + 1), 4);
		break;
	}
}
