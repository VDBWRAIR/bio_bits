*** Settings ***
Library         Process
Library         OperatingSystem
Library         Collections 
Suite Teardown          Terminate All Processes    

*** Variables ***
${in1} =       tests/testinput/780.1.vcf 
${in2} =       tests/testinput/780.2.vcf 
${expected1} =          tests/expected/780.1.out.vcf
${expected2} =          tests/expected/780.2.out.vcf
${vcallin} =    tests/testinput/947.vcf 
${vcallexpected} =      tests/expected/947.vcalls
${existsin} =   tests/testinput/small.vcf
${existsexpected} =      tests/expected/exists.small.vcf

*** Test Cases ***
TestVCFCatDiff
        ${diff_result} =    Run Process     vcfcat  diff    ${in1}  ${in2}  --tag   DP      --threshold     30
        Log To Console      ${diff_result.stderr} 
        Should Be Equal As Integers         ${diff_result.rc}       0
        ${expected_contents} =  Get File        ${expected1}
        Log     ${diff_result.stdout}
        Lists Should Be Equal           ${diff_result.stdout.split()}      ${expected_contents.split()}   
        Log To Console          "Now trying reverse order"
        ${diff_result} =    Run Process     vcfcat  diff    ${in2}  ${in1}  --tag   DP      --threshold     30
        Log To Console      ${diff_result.stderr} 
        Should Be Equal As Integers         ${diff_result.rc}       0
        ${expected_contents} =  Get File        ${expected2}
        Lists Should Be Equal   ${diff_result.stdout.split()}      ${expected_contents.split()}   
        
TestVCFCatVCalls
        ${vcall_result} =    Run Process     vcfcat  vcall    ${vcallin}     --csv
        Log To Console      ${vcall_result.stderr} 
        Should Be Equal As Integers         ${vcall_result.rc}       0 
        ${expected_contents} =  Get File        ${vcallexpected}
        Should Be Equal As Strings     ${vcall_result.stdout.split()}          ${expected_contents.split()} 

TestVCFCatExists
        ${exists_result} =    Run Process     vcfcat  exists    ${existsin}      --tag   ALT
        Log To Console      ${exists_result.stderr} 
        Should Be Equal As Integers         ${exists_result.rc}       0 
        ${expected_contents} =  Get File        ${exists_expected}
        Log     ${exists_result.stdout}
        Should Be Equal As Strings     ${exists_result.stdout.split()}          ${expected_contents.split()} 
        ${exists_result} =    Run Process     vcfcat  exists    ${existsin}      --tag   ALT     --count
        Should Be Equal As Integers     ${exists_result.stdout}  2

