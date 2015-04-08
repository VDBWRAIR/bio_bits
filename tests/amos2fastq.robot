*** Settings ***
Library         Process
Library         OperatingSystem
Suite Teardown          Terminate All Processes    

*** Variables ***
@{EXPECTED}       expected_fooCTG0.fastq           expected_fooCTG1.fastq         expected_fooCTG2.fastq          expected_fooCTG3.fastq          expected_fooCTG4.fastq
@{RESULT}       contig-0.fastq          contig-1.fastq  contig-2.fastq  contig-3.fastq
${in_fastq1} =          foo1.fastq
${in_fastq2} =          foo2.fastq
${in_amos} =            foo.afg

*** Test Cases ***
TestAmos2Fastq
    ${process_result} =         Run Process     amos2fastq      ${in_fastq1}            ${in_fastq2}    ${in_amos}
    # Check system exited  correctly
    Should Be Equal As Integers         ${process_result.rc}        0 
    # Check output
    Should Not Contiain         ${process_result.stdout}        Error

# This can become generic enough to load variables from a file (expected and actual filenames)
# Then assert equality for all of them.

    : FOR ${i}          IN RANGE        0       4
        ${expected} =       Get From List   ${EXPECTED}       ${i}
        ${actual} =     Get From List   ${ACTUAL}       ${I}
        Log     ${expected}
        Log     ${actual}
        File Should Exist       ${actual}
        File Should Not Be Empty        ${actual}
        ${actual_contents} =    Get File        ${actual}
        ${expected_contents} =  Get File        ${expected}
        Should Be Equal As Strings =    ${actual_contents}      ${expected_contents}
        
