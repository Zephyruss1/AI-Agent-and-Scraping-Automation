from streamlit import st


def display_header():
    st.title("Web Search Tool")
    st.subheader("Find Email Addresses and Job Titles")


def input_spreadsheet_link():
    return st.text_input("Enter the link to your spreadsheet:")


def select_model():
    model_choice = st.selectbox(
        "Select the model you want to use:",
        options=["ChatGPT [Browser-Use]", "Perplexity"],
    )
    return model_choice


def select_browser_model():
    browser_model = st.selectbox(
        "Select the browser model:", options=["gpt-4o", "gpt-4o-mini"]
    )
    return browser_model


def select_perplexity_model():
    perplexity_model = st.selectbox(
        "Select the Perplexity model:", options=["sonar-pro", "sonar-deep-research"]
    )
    return perplexity_model


def select_search_purpose():
    purpose_choice = st.selectbox(
        "Select the purpose of your search:",
        options=["Find email addresses", "Find job titles"],
    )
    return purpose_choice


def display_results(results):
    st.subheader("Results")
    for result in results:
        st.write(result)


def display_error(message):
    st.error(message)
