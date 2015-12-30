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
${in_annotation_tab_with_cds} =    tests/testinput/${in_genbank_id}.annotation.withcds.csv

*** Test Cases ***
degen_regions Expected Output Genbank File
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

degen_regions Expected Output Genbank Accession
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

degen_regions Expected Output Tab File
    ${process_result} =    Run Process    degen_regions    -i    ${in_fasta}    -o    ${ACTUAL}    --tab-file    ${in_annotation_tab}    --cds    84,10262

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

degen_regions Expected Output Tab File has cds
    ${process_result} =    Run Process    degen_regions    -i    ${in_fasta}    -o    ${ACTUAL}    --tab-file    ${in_annotation_tab_with_cds}

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

degen_regions command line cds overrides
    ${process_result} =    Run Process    degen_regions    -i    ${in_fasta}    -o    ${ACTUAL}    --gb-id    ${in_genbank_id}    --cds    1,10600

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
    Should Be Not Equal As Strings    ${expected_contents}    ${actual_contents} 
