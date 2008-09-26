/*=========================================================================
Program:   Diffusion Weighted MRI Reconstruction
Module:    $RCSfile: port_ini.h,v $
Language:  C
Author:    $Author: bjian $
Date:      $Date: 2006/07/16 08:15:18 $
Version:   $Revision: 1.3 $
=========================================================================*/

/*
From: Martin Sjodin (msn@ausys.se)
Subject: Reading a Microsoft ini file 
 
  
View this article only 
Newsgroups: comp.unix.programmer
Date: 1995/05/30 
 

In our project we use HP UNIX and programming in C++. The initialisation 
data to the system is given on a diskette with Microsoft ini 
format(text). In Windows there is a command GetPrivateProfileString that 
can be used to read the ini file. The command reads a specific line in a 
specified section and sets a optional variable. My question is: Is there 
a C or C++ program that works in a similar way? There must be! I need the 
source code, and it has to be written in C or C++. Do you have a similar 
program (C or C++ source code) and the corrsponding ini file it should 
help me a lot.

Below there is a short example of a Microsoft ini file. 


Please CC your response to me via email.

Thanks, 
Martin

email: msn@ausys.se


; This is an example of a Microsoft ini file that has to be decoded
; remark

[countries]
France = F
Sweden = S

[France]
Red = 1
black = 2

[Sweden]
Red = 3
black = 4
Message 2 in thread 
From: Paul Mainville (paul@argus.CAM.ORG)
Subject: Re: Reading a Microsoft ini file 
 
  
View this article only 
Newsgroups: comp.unix.programmer
Date: 1995/05/31 
 

I have included two files from an article in Dr Dobbs Journal
from the March 1994 issue. 

Modifications from the article: use same names as for MS-Win API and
added a WritePrivateProfileInt().

Paul

============================================================================
Paul Mainville                         email: paul@argus.cam.org 
Boisbriand, Quebec, Canada
============================================================================


*/

/*F    func name: port_ini.h
 *
 *         descr: Portable Routines for writing PRIVATE PROFILE STRINGS
 *                From DDJ March 1994.
 *
 *      designed: See Article
 * last modified: Paul Mainville 94/07/01
 *
 * $Log: port_ini.h,v $
 * Revision 1.3  2006/07/16 08:15:18  bjian
 * *** empty log message ***
 *
 * Revision 1.2  2006/06/20 16:20:09  bjian
 * *** empty log message ***
 *
 * Revision 1.1  2006/03/23 05:27:44  bjian
 * *** empty log message ***
 *
 * Revision 1.3  2005/03/17 13:22:49  bjian
 * *** empty log message ***
 *
 * Revision 1.2  2004/02/09 20:57:08  bjian
 * *** empty log message ***
 *
 * Revision 1.1  2004/02/09 20:33:21  bjian
 * *** empty log message ***
 *
 */

#ifndef _PORT_INI_H_
#define _PORT_INI_H_


#include <string.h>
#include <stdio.h>


# ifdef __cplusplus
  extern "C"{
# endif

#define MAX_LINE_LENGTH    180

int GetPrivateProfileInt(const char *, const char *, int, const char *);
int GetPrivateProfileString(const char *, const char *, const char *, char *,int, const char *);
int WritePrivateProfileString(char *, char *, char *, char *);
int WritePrivateProfileInt(char *, char *, int, char *);


# ifdef __cplusplus
  }
# endif

#endif
/* End of port_ini.h */
