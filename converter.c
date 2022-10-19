/****************************************************************/
/*	convbin2asc		01.08.2004			*/
/****************************************************************/
/*	Short Description :                                     */
/*	AU program to convert 1r and/or 1i data file into an	*/
/*	ascii table containing point number, Hz, ppm and	*/
/*	intensity values.					*/
/****************************************************************/
/*	Keywords :                                              */
/*	ascii conversion					*/
/****************************************************************/
/*	Description/Usage :                                     */
/*	AU program to convert 1r and/or 1i data file into an	*/
/*	ascii table containing point number, Hz, ppm and	*/
/*	intensity values. To get the imaginary data points as	*/
/*	as well, start the AU program with: "convbin2asc i"	*/
/*	The output file is stored in the same directory as the	*/
/*	spectrum. It can be used for calculations in spread-	*/
/*	sheet programs or other third-party software.		*/
/****************************************************************/
/*	Author(s) :                                             */
/*	Name		: Clemens Anklin			*/
/*	Organisation	: Bruker BioSpin			*/
/*	Email		: cga@bruker.com			*/
/****************************************************************/
/*	Name            Date    Modification:			*/
/*	cga		040801  created				*/
/*	eng		050201  byte order check added		*/
/****************************************************************/
/*
$Id: convbin2asc,v 1.1.2.2 2005/02/01 14:55:10 eng Exp $
*/

#include <fcntl.h>

char infile[PATH_MAX], outfile[PATH_MAX], also1i[2], *inbuf1;
int *inbufint1, pparmode, si, i, fpin, intens, bytordp;
unsigned bytetoread;
float offset, hzppt;
double sf, swp, hz, ppm,  hzoffset, ppmppt;
FILE *fpout;

GETCURDATA

FETCHPARS("PPARMOD",&pparmode);
if (pparmode != 0)
{
  Proc_err(DEF_ERR_OPT, "convbin2asc only works for 1D data.");
  ABORT
}

also1i[0] = '\0';
sscanf(cmd,"%s",also1i);

FETCHPARS("BYTORDP",&bytordp);
FETCHPARS("SI",&si);
FETCHPARS("OFFSET",&offset);
FETCHPARS("SF",&sf);
FETCHPARS("SW_p",&swp);
FETCHPARS("HzpPT",&hzppt);

ppmppt=(double)hzppt/sf;
hzoffset=(double)offset*sf;
bytetoread = si*sizeof(int);

/*
Proc_err(DEF_ERR_OPT, "offset %.4f, sf %.7f, SW_p %.3f,\n"
"HZpPT %.2f, ppmppt %.4f, hzoffset %.2f, bytetoread %d",
offset, sf, swp, hzppt, ppmppt, hzoffset, bytetoread );
*/

(void)sprintf(infile,"%s/data/%s/nmr/%s/%d/pdata/%d/1r",
disk,user,name,expno,procno);

(void)sprintf(outfile,"%s/data/%s/nmr/%s/%d/pdata/%d/ascii-spec.txt",
disk,user,name,expno,procno);

if ( (fpin=open(infile,O_RDONLY)) == (-1) )
{
  Perror(DEF_ERR_OPT,infile);
  ABORT
}
if ( (fpout=fopen(outfile,"wb")) == NULL)
{
  Perror(DEF_ERR_OPT,outfile);
  ABORT;
}

if ( (inbuf1 = calloc (bytetoread,sizeof(char*))) == NULL)
{
  Perror(DEF_ERR_OPT,infile);
  close (fpin);
  fclose (fpout);
  ABORT
}
if ( (i = read (fpin,inbuf1,bytetoread)) != bytetoread )
{
  Perror(DEF_ERR_OPT,infile);
  fclose (fpout);
  close (fpin);
  ABORT
}
close(fpin);
local_swap4(inbuf1,bytetoread,bytordp);

inbufint1 = (int *)inbuf1;

for (i=0; i < si; i++)
  {
  intens=inbufint1[i];
  ppm=offset-ppmppt*(i+0.5);
  hz=hzoffset-hzppt*(i+0.5);
/*    #fprintf(fpout,"%i, %d , %.3f, %.4f\n",i+1,intens,hz,ppm);
*/
    fprintf(fpout,"%.4f	%d\n",ppm, intens);
  }

if (also1i[0] != '\0')
{
  (void)sprintf(infile,"%s/data/%s/nmr/%s/%d/pdata/%d/1i",
disk,user,name,expno,procno);

  if ( (fpin=open(infile,O_RDONLY)) == (-1) )
  {
    Perror(DEF_ERR_OPT,infile);
    ABORT
  }
  if ( (i = read (fpin,inbuf1,bytetoread)) != bytetoread )
  {
    Perror(DEF_ERR_OPT,infile);
    fclose (fpout);
    close (fpin);
    ABORT
  }
  close(fpin);
  local_swap4(inbuf1,bytetoread,bytordp);
  for (i=0; i < si; i++)
  {
    intens=inbufint1[i];
    ppm=offset-ppmppt*(i+0.5);
    hz=hzoffset-hzppt*(i+0.5);
/*    #fprintf(fpout,"%i, %d , %.3f, %.4f\n",i+1,intens,hz,ppm);
*/
    fprintf(fpout,"%.4f	%d\n",ppm, intens);
  }
}

fclose(fpout);

QUIT