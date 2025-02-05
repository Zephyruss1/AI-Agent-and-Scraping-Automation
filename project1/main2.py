from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from webdriver_manager.chrome import ChromeDriverManager
import time
import pyautogui
import subprocess


def connect_to_protonvpn():
    ovpn_config = "C:\Users\ekber\OpenVPN\config\.ovpn-main\Vietnam\IPS_118.71.213.89_udp_1759.ovpn"
    auth_file = "/path/to/auth.txt"

    try:
        subprocess.run(["sudo", "openvpn", "--config", ovpn_config, "--auth-user-pass", auth_file], check=True)
    except subprocess.CalledProcessError as e:
        print(f"OpenVPN connection failed: {e}")

connect_to_protonvpn()


def run_command():
    # URL to search for
    url_to_search = "https://www.marsbahisonline.com"

    # Set up Chrome options
    options = Options()
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    options.add_argument("--ignore-certificate-errors")
    options.add_argument("--headless")  # Uncomment to run in headless mode

    # Initialize the WebDriver
    driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()), options=options)
    while True:
        try:
            # Step 1: Open DuckDuckGo
            driver.get("https://www.duckduckgo.com")
            time.sleep(2)  # Wait for the page to load

            # Step 2: Find the search box and enter the URL
            search_box = driver.find_element(By.NAME, "q")
            search_box.send_keys(url_to_search)
            search_box.send_keys(Keys.RETURN)  # Press Enter to search
            time.sleep(2)  # Wait for the search results to load

            # Step 3: Click on the first search result
            # Wait for the first search result to be present
            first_result = WebDriverWait(driver, 10).until(
                EC.presence_of_element_located((By.CSS_SELECTOR, "h2 a"))  # Correct selector for DuckDuckGo
            )
            first_result.click()  # Click on the first result
            time.sleep(5)  # Wait for the page to load

            print("Successfully navigated to the URL via DuckDuckGo search.")

        except Exception as e:
            print(f"An error occurred: {e}")

        finally:
            # Close the browser
            import time
            time.sleep(10)
            driver.quit()