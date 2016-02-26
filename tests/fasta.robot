*** Settings ***
Library             Process
Library             OperatingSystem
Library             Collections 
Library             String
Suite Teardown      Terminate All Processes    

*** Variables ***
${in_fasta} =                       ${CURDIR}/testinput/col.fasta
${actual} =                         ${CURDIR}/out.fasta
${expected} =                       ${CURDIR}/expected/singleline.fasta
${test_directory} =                 ${CURDIR}/output

*** Test Cases ***
fasta reads from stdin
    ${process_result} =             Run Process                     cat ${in_fasta} | fasta - > ${actual}   shell=True

    # Check system exited    correctly
    Should Be Equal As Integers     ${process_result.rc}            0 

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

    # Check output
    Should Not Contain              ${process_result.stdout}        Error
    File Should Exist               ${actual}
    File Should Not Be Empty        ${actual}
    ${actual_contents} =            Get File                        ${actual}
    ${expected_contents} =          Get File                        ${expected}
    Should Be Equal As Strings      ${expected_contents}            ${actual_contents} 
fasta used in shell pipeline
    ${process_result} =             Run Process                     fasta ${in_fasta} | grep -v '>' | grep -Eo '[Aa]' | wc -l    shell=True

    # Check system exited    correctly
    Should Be Equal As Integers     ${process_result.rc}            0 

    # Check output
    Should Be Equal                 ${process_result.stdout}        160
fasta wraps sequences
    ${process_result} =             Run Process                     fasta ${in_fasta} | fasta --wrap - > ${actual}    shell=True

    # Check system exited    correctly
    Should Be Equal As Integers     ${process_result.rc}            0 

    # Check output
    ${actual_contents} =            Get File                        ${actual}
    ${expected_contents} =          Get File                        ${in_fasta}
    Should Be Equal As Strings      ${expected_contents}            ${actual_contents} 
fasta splits all identifiers to new files
    Create Directory                ${test_directory}
    Empty Directory                 ${test_directory}
    ${process_result} =             Run Process                     fasta --split ${in_fasta}    shell=True    cwd=${test_directory}

    # Check system exited    correctly
    Should Be Equal As Integers     ${process_result.rc}            0 
    @{idents} =                     Create List                     sequence1                       sequence2____________________________
    :FOR    ${id}    in    @{idents}
    \       ${contents} =    Get File    ${test_directory}/${id}.fasta
    \       ${numlines} =    Get Line Count  ${contents}
    \       Should Be Equal As Integers    2    ${numlines}
fasta splits works with stdin
    Create Directory                ${test_directory}
    Empty Directory                 ${test_directory}
    ${process_result} =             Run Process                     cat ${in_fasta} | fasta --split -    shell=True    cwd=${test_directory}

    # Check system exited    correctly
    Should Be Equal As Integers     ${process_result.rc}            0 
    @{idents} =                     Create List                     sequence1                       sequence2____________________________
    :FOR    ${id}    in    @{idents}
    \       ${contents} =    Get File    ${test_directory}/${id}.fasta
    \       ${numlines} =    Get Line Count  ${contents}
    \       Should Be Equal As Integers    2    ${numlines}
