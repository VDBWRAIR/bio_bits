*** Settings ***
Library         Process
Library         OperatingSystem
Suite Teardown          Terminate All Processes    

*** Keywords ***
Should Be Equal As Files        [Arguments]     ${file1}        ${file2}
        ${contents1} =          Get File        ${file1}
        ${contents2} =          Get File        ${file2}
        Log To Console          ${contents1}
        Log To Console          ${contents2}
        Should Be Equal as Strings      ${contents1}    ${contents2}

*** Variables ***
${in1} =        tests/testinput/out.samtext
${out1} =       chr1.group.fq
${out2} =       chr2.group.fq

*** Test Cases ***
TestParseRefs
        ${result} =     Run Process     parse_contigs      ${in1}
        Log To Console      ${result.stderr} 
        Should Be Equal As Integers         ${result.rc}       0
        Should Be Equal As Files        tests/expected/${out1}          ${out1}
        Should Be Equal As Files        tests/expected/${out2}          ${out2}

