*** Settings ***
Library    Process
Library    OperatingSystem
Library    Collections 
Suite Teardown    Terminate All Processes    

*** Variables ***
${EXPECTED} =    tests/testinput/ctl_expected.tsv
${ACTUAL} =    robotout.tsv
${in_fasta} =    tests/Den4_MAAPS_TestData16.fasta 
${in_genbank} =    tests/testinput/sequence.gb
${in_genbank_id} =    KJ189367
${in_annotation_tab} =    tests/testinput/${in_genbank_id}.annotation.csv

*** Test Cases ***
Expected Output Genbank File
    ${process_result} =    Run Process    degen_regions    -i    ${in_fasta}    -o    ${ACTUAL}    --gb-file    ${in_genbank}

    # Check system exited    correctly
    Should Be Equal As Integers    ${process_result.rc}    0 
    Log To Console    ${process_result.stdout}
    Log To Console    ${process_result.stderr}

    # Check output
    Should Not Contain    ${process_result.stdout}    Error
    File Should Exist    ${actual}
    File Should Not Be Empty    ${actual}
    ${actual_contents} =    Get File    ${actual}
    ${expected_contents} =    Get File    ${expected}
    Should Be Equal As Strings    ${expected_contents}    ${actual_contents} 

Expected Output Genbank Accession
    ${process_result} =    Run Process    degen_regions    -i    ${in_fasta}    -o    ${ACTUAL}    --gb-id    ${in_genbank_id}

    # Check system exited    correctly
    Should Be Equal As Integers    ${process_result.rc}    0 
    Log To Console    ${process_result.stdout}
    Log To Console    ${process_result.stderr}

    # Check output
    Should Not Contain    ${process_result.stdout}    Error
    File Should Exist    ${actual}
    File Should Not Be Empty    ${actual}
    ${actual_contents} =    Get File    ${actual}
    ${expected_contents} =    Get File    ${expected}
    Should Be Equal As Strings    ${expected_contents}    ${actual_contents} 

Expected Output Tab File
    ${process_result} =    Run Process    degen_regions    -i    ${in_fasta}    -o    ${ACTUAL}    --tab-file    ${in_annotation_tab}

    # Check system exited    correctly
    Should Be Equal As Integers    ${process_result.rc}    0 
    Log To Console    ${process_result.stdout}
    Log To Console    ${process_result.stderr}

    # Check output
    Should Not Contain    ${process_result.stdout}    Error
    File Should Exist    ${actual}
    File Should Not Be Empty    ${actual}
    ${actual_contents} =    Get File    ${actual}
    ${expected_contents} =    Get File    ${expected}
    Should Be Equal As Strings    ${expected_contents}    ${actual_contents} 
