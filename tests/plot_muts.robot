*** Settings ***
Library             Process
Library             OperatingSystem
Library             Collections 
Library             String
Suite Teardown      Terminate All Processes    


*** Variables ***
${png} =         tests/muts.png
${csv} =         ${png}.csv
${html} =        ${png}.html
${expected} =   tests/expected/out.ha.png.csv


*** Test Cases ***
plot_muts test
	${process_result} =             Run Process         plot_muts --query tests/testinput/ha/query.ha.fasta --refs tests/testinput/ha/refall.ha.fasta --out ${png}   shell=True
        Log To Console       ${process_result.stdout}
        Log To Console        ${process_result.stderr}
        Should Be Equal As Integers         ${process_result.rc}        0 
        File Should Exist       ${png}
        File Should Not Be Empty        ${png}
        File Should Exist       ${csv}
        File Should Not Be Empty        ${csv}
        ${actual_contents} =            Get File                        ${csv}
        ${expected_contents} =          Get File                        ${expected} 
        Should Be Equal As Strings      ${expected_contents}            ${actual_contents} 

plot_muts html 
	${process_result} =             Run Process         plot_muts --query tests/testinput/ha/query.ha.fasta --refs tests/testinput/ha/refall.ha.fasta --html --out ${png}   shell=True
        Log To Console       ${process_result.stdout}
        Log To Console        ${process_result.stderr}
        Should Be Equal As Integers         ${process_result.rc}        0 
        File Should Exist       ${png}
        File Should Not Be Empty        ${png}
        File Should Exist       ${csv}
        File Should Not Be Empty        ${csv}
        File Should Exist       ${html}
        File Should Not Be Empty        ${html}

