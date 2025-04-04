C     Originallly Coded:  D. Gilmore
C           Modified By:  D. Gilmore ---> 19OCT95
C     -----------------------------------------------------
C     L. ADAMOWICZ RESEARCH / UNIVERSITY OF ARIZONA / CHEM.
C     -----------------------------------------------------
      SUBROUTINE BIOUT(STRING)
C
C     BIOUT prints string output to both the OUTPUT.TXT
C         file and the screen
C
      CHARACTER*(*) STRING
C
      WRITE(6,*) STRING
      WRITE(7,*) STRING
C
      RETURN
      END
