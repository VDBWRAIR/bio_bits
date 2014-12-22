from subprocess import Popen

def mpileup(self, *args):
    '''
    Call samtools mpileup with args
    Accepts same args as mpileup and returns iterater of returned mpileup rows
    '''
