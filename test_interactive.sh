#!/bin/bash

# Test the quick_start.sh interactively
echo "Testing quick_start.sh with manual input simulation..."

# Create a simple test
echo "1. First test - checking help (option 4)"
echo "4" | timeout 10 bash quick_start.sh

echo ""
echo "2. Test completed. If you want to test manually, run:"
echo "   bash quick_start.sh"
echo "   and choose option 4 to see the help"