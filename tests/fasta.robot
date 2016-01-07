*** Settings ***
Library             Process
Library             OperatingSystem
Library             Collections 
Library             String
Suite Teardown      Terminate All Processes    

*** Variables ***
${in_fasta} =                       tests/testinput/col.fasta
${actual} =                      tests/out.fasta
${expected} =                 tests/expected/singleline.fasta

*** Test Cases ***
fasta reads from stdin
    ${process_result} =             Run Process                     cat ${in_fasta} | fasta - > ${actual}   shell=True

    # Check system exited    correctly
    Should Be Equal As Integers     ${process_result.rc}            0 
    Log To Console                  ${process_result.stdout}
    Log To Console                  ${process_result.stderr}

    # Check output
    Should Not Contain              ${process_result.stdout}        Error
    File Should Exist               ${actual}
    File Should Not Be Empty        ${actual}
    ${actual_contents} =            Get File                        ${actual}
    ${expected_contents} =          Get File                        ${expected}
    Should Be Equal As Strings      ${expected_contents}            ${actual_contents} 
fasta reads from file
    ${process_result} =             Run Process                     fasta ${in_fasta} > ${actual}           shell=True

    # Check system exited    correctly
    Should Be Equal As Integers     ${process_result.rc}            0 
    Log To Console                  ${process_result.stdout}
    Log To Console                  ${process_result.stderr}

    # Check output
    Should Not Contain              ${process_result.stdout}        Error
    File Should Exist               ${actual}
    File Should Not Be Empty        ${actual}
    ${actual_contents} =            Get File                        ${actual}
    ${expected_contents} =          Get File                        ${expected}
    Should Be Equal As Strings      ${expected_contents}            ${actual_contents} 
