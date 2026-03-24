# Data Model Classes

class User:
    def __init__(self, user_id, username, email):
        self.user_id = user_id
        self.username = username
        self.email = email

class Product:
    def __init__(self, product_id, name, price):
        self.product_id = product_id
        self.name = name
        self.price = price

class Order:
    def __init__(self, order_id, user, product, quantity):
        self.order_id = order_id
        self.user = user
        self.product = product
        self.quantity = quantity