/*=========================================================================
Program:   Diffusion Weighted MRI Reconstruction
Module:    $RCSfile: port_ini.c,v $
Language:  C
Author:    $Author: bjian $
Date:      $Date: 2006/07/16 08:15:18 $
Version:   $Revision: 1.3 $
=========================================================================*/


/*F    func name: port_ini.c
 *
 *         descr: Portable Routines for writing PRIVATE PROFILE STRINGS
 *                From DDJ March 1994.
 *
 *      designed: See Article
 * last modified: Paul Mainville 94/07/01
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "port_ini.h"

#pragma warn -wpia      /* BC++ Ignore possible incorrect assignment */

/*F    func name: read_line()
 *
 *         descr: Read a line from the file.
 *
 *    parameters: FILE *   fp    pointer to the file to read from.
 *                char *   bp    pointer to the copy buffer.
 *
 *       returns: int      TRUE of success FALSE otherwise
 */

int read_line(FILE *fp, char *bp)
{
   char  c = '\0';
   int   i = 0;

   while ((c = getc(fp)) != '\n' && c!='\r')
   {
      if (c == EOF)
         return 0;

      bp[i] = c;
      i++;
   }

   bp[i] = '\0';
   //printf("[read_line] %s\n", bp);
   return 1;
}

/*F    func name: GetPrivateProfileInt()
 *
 *         descr: Get integer from INI file.
 *
 *    parameters: char *section     name of the section to search for
 *                char *entry       name of the entry to find the value of
 *                int   def         default value if not found
 *                char *file_name   name of the .ini file to read from
 *
 *       returns: int               the value read.
 */

int GetPrivateProfileInt(const char *section, const char *entry, int def, const char *file_name)
{
   FILE  *fp = fopen(file_name, "r");
   char  buff[MAX_LINE_LENGTH];
   char  *ep;
   char  t_section[MAX_LINE_LENGTH];
   char  value[6];
   int   len = strlen(entry);
   int   i;

   //!  printf("%s %s %s\n", section, entry, file_name);

   if (!fp)
      return -1;

   //! printf("%s %s %s\n", section, entry, file_name);

   sprintf(t_section, "[%s]", section);

   do       /* Search until section is found or EOF */
   {
      if (!read_line(fp, buff))
      {
         fclose(fp);
         return(def);
      }
   } while (strcmp(buff, t_section));

   do       /* Section found - search entry - stop search at end of section */
   {
      if (!read_line(fp, buff) || (buff[0] == '\n'))
      {
         fclose(fp);
         return(def);
      }
   } while (strncmp(buff, entry, len));


   /*  ep = strrchr(buff, '=');            /* parse out '=' sign */
   /* ep++; */
   //! printf("0--%s\n",buff);   
   ep = strtok(buff, "\t =\n\r");
   //! printf("1--%s\n",ep);
   ep = strtok(NULL, "\t =\n\r");
   //! printf("2--%s\n",ep);

   if (!strlen(ep))
      return(def);

   for (i = 0; isdigit(ep[i]); i++)    /* copy only numbers */
      value[i] = ep[i];
   
   value[i] = '\0';

   //! printf("v--%s\n",value);
   fclose(fp);
   return(atoi(value));
}

/*F    func name: GetPrivateProfileString()
 *
 *         descr: Get string from INI file.
 *
 *    parameters: char *section     name of the section to search for
 *                char *entry       name of the entry to find the value of
 *                char *def         default value if not found
 *                char *buffer      pointer to the buffer to write into
 *                int   buffer_len  max number of char to copy
 *                char *file_name   name of the .ini file to read from
 *
 *       returns: int               number of char copied
 */

int GetPrivateProfileString(const char *section, const char *entry, const char *def, 
                            char *buffer, int buffer_len, const char *file_name)
{

   FILE  *fp; 
   char  buff[MAX_LINE_LENGTH];
   char  *ep;
   char  t_section[MAX_LINE_LENGTH];
   int   len = strlen(entry);
   int   sec_len;
   int flag;

   fp  = fopen(file_name, "r");


   //printf("getstr0 %s %s %s\n", section, entry, file_name);

   if (!fp)
      return 0;

   // printf("getstr1 %s %s %s\n", section, entry, file_name);

   sprintf(t_section, "[%s]", section);
   sec_len = strlen(t_section);
   t_section[sec_len] = '\0';

   //printf("getstr1.1 %s %s %s\n", t_section, entry, file_name);


   do       /* Search until section is found or EOF */
   {
     //printf("[do..while..loop] %s %s %s\n", t_section, entry, file_name);
      
     // printf("???? buff = %s \n", buff);
      flag = read_line(fp,buff);
      //printf("before flag: buff = %s\n",buff);
      //printf("ead [%s]\n", buff); 
      if (!flag) 
      //if (!read_line(fp, buff))
      {
	//printf("getstr1.5 %s %s\n", buff, t_section);
         fclose(fp);
         if (def) strncpy(buffer, def, buffer_len);
         return(strlen(buffer));
      }
      //printf("read [%s]\n",buff);
      //printf("try to compare [%s]\n",t_section);
   } while (strcmp(buff, t_section));

   //! printf("getstr2 %s %s %s\n", section, entry, file_name);

   do       /* Section found - search entry - stop search at end of section */
   {
      if (!read_line(fp, buff) || (buff[0] == '\n'))
      {
         fclose(fp);
         if (def) strncpy(buffer, def, buffer_len);
         return(strlen(buffer));
      }
   } while (strncmp(buff, entry, len));

   //! printf("getstr3 %s %s %s\n", section, entry, file_name);

   /*  ep = strrchr(buff, '=');            /* parse out '=' sign */
   /* ep++; */
   
   ep = strtok(buff, "\t =\n\r");

   //! printf("getstr4 %s %s %s ep %s buffer %s \n", section, entry, file_name, ep, buffer);

   ep = strtok(NULL, "\t =\n\r");

   //! printf("getstr5 %s %s %s %p ep %s buffer %s %p \n", section, entry, file_name, file_name, ep, buffer, buffer);

   if (!ep || strlen(ep) == 0)
     {
       fclose(fp);
       //!printf("getstr5 %s %s %s %p ep %s buffer %s %p \n", section, entry, file_name, file_name, ep, buffer, buffer);
       if (def) strncpy(buffer, def, buffer_len);
       return(strlen(buffer));
     }   
   strncpy(buffer, ep, strlen(ep));
   //! printf("getstr6 %s %s %s %p ep %s buffer %s %p\n", section, entry, file_name, file_name, ep, buffer, buffer);

   buffer[strlen(ep)] = '\0';
   //! printf("getstr7 %s %s %s %p ep %s buffer %s %p\n", section, entry, file_name, file_name, ep, buffer, buffer);


   fclose(fp);
   return(strlen(buffer));
}


/*F    func name: WritePrivateProfileString()
 *
 *         descr: Write string into INI file
 *
 *    parameters: char *section     name of the section to search for
 *                char *entry       name of the entry to find the value of
 *                char *buffer      pointer to the buffer to write into
 *                char *file_name   name of the .ini file to read from
 *
 *       returns: int               TRUE if succesful, FALSE otherwise
 */

int WritePrivateProfileString(char *section, char *entry, char *buffer,
                              char *file_name)
{
   FILE  *rfp, *wfp;
   char  tmp_name[15];
   char  buff[MAX_LINE_LENGTH];
   char  t_section[MAX_LINE_LENGTH];
   int   len = strlen(entry);

   mkstemp(tmp_name);
   sprintf(t_section, "[%s]", section);

   if (!(rfp = fopen(file_name, "r")))    /* If the INI file does not exist */
   {                                      /* Then make one */
      if (!(wfp = fopen(file_name, "w")))
         return 0;

      fprintf(wfp, "%s\n", t_section);
      fprintf(wfp, "%s=%s\n", entry, buffer);
      fclose(wfp);
      return 1;
   }

   if (!(wfp = fopen(tmp_name, "w")))
   {
      fclose(wfp);
      return 0;
   }

   /* Move through the file one line at a time until a section is matched */
   /* or until EOF. Copy to temp file as it is read. */

   do
   {
      if (!read_line(rfp, buff))
      {
         /* Failed to find section. Add one to the end. */

         fprintf(wfp, "\n%s\n", t_section);
         fprintf(wfp, "%s=%s\n", entry, buffer);

         fclose(rfp);
         fclose(wfp);
         unlink(file_name);
         rename(tmp_name, file_name);
         return 1;
      }

      fprintf(wfp, "%s\n", buff);
   } while (strcmp(buff, t_section));

   /* Now that the section has been found, find the entry. Stop searching */
   /* upon leaving the section's area. Copy the file as it is read, and   */
   /* create an entry if one is not found.                                */

   while(1)
   {
      if (!read_line(rfp, buff))
      {
         /* EOF without an entry so make one */
         fprintf(wfp, "%s=%s\n", entry, buffer);

         fclose(rfp);
         fclose(wfp);
         unlink(file_name);
         rename(tmp_name, file_name);
         return 1;
      }

      if (!strncmp(buff, entry, len) || buff[0] == '\0')
         break;

      fprintf(wfp, "%s\n", buff);
   }

   if (buff[0] == '\0')
   {
      fprintf(wfp, "%s=%s\n", entry, buffer);

      do
      {
         fprintf(wfp, "%s\n", buff);
      } while (read_line(rfp, buff));
   }
   else
   {
      fprintf(wfp, "%s=%s\n", entry, buffer);

      while (read_line(rfp, buff))
         fprintf(wfp, "%s\n", buff);
   }
      
   fclose(rfp);
   fclose(wfp);
   unlink(file_name);
   rename(tmp_name, file_name);
   return 1;
}
/*F    func name: WritePrivateProfileInt()
 *
 *         descr: Write integer into INI file
 *
 *    parameters: char *section     name of the section to search for
 *                char *entry       name of the entry to find the value of
 *                int   value       value to write
 *                char *file_name   name of the .ini file to read from
 *
 *       returns: int               TRUE if succesful, FALSE otherwise
 */

int WritePrivateProfileInt(char *section, char *entry, int value,
                              char *file_name)
{
   char  buff[10];

   /* itoa(value, buff, 10); itoa nonstandard */
   sprintf(buff, "%d", value);

   return WritePrivateProfileString(section, entry, buff, file_name);
}

/* End of port_ini.c */
