from pyautogui import press
import time


'''
script that toggles the numlock key every 9 mins (540secs) 
to keep screen from locking
'''
start_time = time.time()

print("Keeping the Screen Unlocked.")



while True:
		press('numlock')
		time.sleep(0.5)
		press('numlock')
		time.sleep(540)
		print("Time running: %s mins" % round((time.time()-start_time)/60,2), end='\r')