from PyPDF2 import PdfMerger, PdfReader
import os

save_file = "LEAP_03_25.pdf"

mergeList = ["Receipt_10-Mar-2025.pdf", "Receipt_17-Mar-2025.pdf", "Receipt_31-Mar-2025.pdf"]
mergeFolder = (os.environ['USERPROFILE']+r"\Downloads")

myMerger = PdfMerger()

for filename in mergeList:
    myMerger.append(PdfReader(os.path.join(mergeFolder, filename)), "rb")

myMerger.write(os.path.join(mergeFolder, save_file))