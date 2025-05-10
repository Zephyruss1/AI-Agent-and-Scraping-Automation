# Streamlit Web Search UI

This project is a Streamlit application designed to facilitate web searches for email addresses and job titles based on author names. It integrates functionality from the `test_for_websearch.py` script to provide a user-friendly interface for users to input their search criteria and view results.

## Project Structure

```
streamlit-websearch-ui
├── src
│   ├── app.py                # Main entry point for the Streamlit application
│   ├── components
│   │   ├── __init__.py       # Marks the components directory as a package
│   │   └── ui_elements.py     # Defines UI elements for the application
│   ├── utils
│   │   ├── __init__.py       # Marks the utils directory as a package
│   │   └── websearch_utils.py  # Utility functions for web search functionality
│   └── assets
│       └── styles.css        # Custom CSS styles for the application
├── requirements.txt           # Lists project dependencies
├── .streamlit
│   └── config.toml           # Configuration settings for the Streamlit application
└── README.md                  # Documentation for the project
```

## Setup Instructions

1. **Clone the Repository**
   ```bash
   git clone <repository-url>
   cd streamlit-websearch-ui
   ```

2. **Install Dependencies**
   It is recommended to use a virtual environment. You can create one using `venv` or `conda`.
   ```bash
   pip install -r requirements.txt
   ```

3. **Run the Application**
   Start the Streamlit application by running:
   ```bash
   streamlit run src/app.py
   ```

## Usage

- Upon launching the application, users will be prompted to enter a link to a spreadsheet containing author names and keywords.
- Users can select the model they wish to use for the search (ChatGPT or Perplexity).
- The application allows users to choose whether they want to find email addresses or job titles.
- Results will be displayed on the interface, showing the found email addresses or job titles.

## Contributing

Contributions are welcome! Please feel free to submit a pull request or open an issue for any enhancements or bug fixes.

## License

This project is licensed under the MIT License. See the LICENSE file for more details.
