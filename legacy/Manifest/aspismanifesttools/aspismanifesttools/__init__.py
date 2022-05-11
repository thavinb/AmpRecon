
import sys
import os

from aspismanifesttools.sanger2manifest import *
from aspismanifesttools.manifest2barcode import *
from aspismanifesttools.manifest2manifest import *
from aspismanifesttools.manifest2taglist import *

#------------------------------------------------- 
def execCmd():
 
 toolName = os.path.basename(sys.argv[0])
 
 supportedCommands = ["sanger2manifest","manifest2barcode","manifest2manifest","manifest2taglist"]
 supportedCalls    = [sanger2ManifestCmd,manifest2barcodeCmd,manifest2manifestCmd,manifest2taglistCmd]

 if len(sys.argv) <= 1:
  printHelp(toolName)
  
 else:
  args = [toolName]
  
  if len(sys.argv) >= 2:
   args.extend(sys.argv[2:])
  
  found = False
  for i in range(0,len(supportedCommands)):
   if sys.argv[1] == supportedCommands[i]:
    args[0] += " {}".format(supportedCommands[i])
    supportedCalls[i](args)
    found = True
    break
    

  if not found:
   print("Unsupported command!")
   printHelp(toolName)

  
#------------------------------------------------- 
def printHelp(toolName):

 print("\n")
 print(toolName)
 print(len(toolName) * "=")
 print("ASPIS Manifest tools")
 print("usage: {} [cmd] [args]\n".format(toolName))
 print("Tools")
 print("-----")
 print("sanger2manifest - creates an ASPIS manifest based on a Sanger run id")
 print("manifest2barcode - uses an ASPIS manifest to generate a barcode file")
 print("manifest2manifest - uses an ASPIS manifest to generate a second phase manifest")
 print("manifest2taglist - uses an ASPIS manifest to a taglist file")
