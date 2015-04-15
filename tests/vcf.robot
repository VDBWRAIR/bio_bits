
*** Settings ***
Library         Process
Library         OperatingSystem
Library         Collections 
Suite Teardown          Terminate All Processes    


*** Variables ***
${in1} =        
${in2} =        
${expected1} =   
${vcallin} =     




*** Test Cases ***
PassForNow
        Log To Console          "functionality missing"
#TestVCFCat 
#        ${diff_result} =    Run Process     vcfcat  diff    ${in1}  ${in2}  --tag   DP      -t     30
#        Log To Console      ${diff_result.stderr} 
#        Should Be Equal As Integers         ${diff_result.rc}       0
#        ${expetced_contents} =  Get File        ${expected1}
#        Should Be Equal As Strings      ${diff_result.stdout}      ${expected_contents}   
#        
#TestVCFCatVCalls
#        ${vcall_result} =    Run Process     vcfcat  vcall    ${vcallin}  -c
#        Log To Console      ${vcall_result.stderr} 
#        Should Be Equal As Integers         ${diff_result.rc}       0
#        Should Be Equal As Integers     ${vcall_result.stdout}          3 
#
#TestVCFCatFilter
#
#TestVCFCatAmbiguous
