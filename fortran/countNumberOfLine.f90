MODULE mod_countNbOfLine
 IMPLICIT NONE
 CONTAINS
    INTEGER FUNCTION countNbOfLine(filename,NbOfLetterInFilename)
        INTEGER                         :: NbOfLetterInFilename
        CHARACTER(NbOfLetterInFilename) :: filename
        CALL EXECUTE_COMMAND_LINE('wc -l <'//filename//'> wc.txt' ) 
        OPEN(unit=1,file='wc.txt') 
        READ(1,*) countNbOfLine
        CLOSE(1)
        CALL EXECUTE_COMMAND_LINE('rm wc.txt')
    END FUNCTION countNbOfLine
END MODULE mod_countNbOfLine
