def extract_proxy_format(proxy_list: list) -> dict:
    """
    Extracts proxies from a list and returns them in a dictionary format.

    Args:
        proxy_list (list): List of proxies in the format 'ip:port:username:password'.

    Returns:
        dict: Dictionary with keys 'ip', 'port', 'username', 'password' and their corresponding values.
    """
    if not proxy_list:
        return []

    formatted_proxy_list = []
    for proxy in proxy_list:
        parts = proxy.split(":")
        if len(parts) == 4:
            proxy_info = {
                "server": f"http://{parts[0]}:{parts[1]}",
                "username": parts[2],
                "password": parts[3],
            }
            formatted_proxy_list.append(proxy_info)
    return formatted_proxy_list
